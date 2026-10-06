{smcl}
{* *! version 1.6.0  26aug2026}{...}
{vieweralsosee "[R] quaids" "help quaids"}{...}
{viewerjumpto "Syntax" "pyquaidsce##syntax"}{...}
{viewerjumpto "Description" "pyquaidsce##description"}{...}
{viewerjumpto "Options" "pyquaidsce##options"}{...}
{viewerjumpto "Examples" "pyquaidsce##examples"}{...}
{viewerjumpto "Stored results" "pyquaidsce##results"}{...}
{viewerjumpto "Author" "pyquaidsce##author"}{...}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{bf:pyquaidsce} {hline 2}}Fast Censored Quadratic Almost Ideal Demand System (QUAIDS) Estimation in Stata via Python{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 18 2}
{cmd:pyquaidsce} {it:sharelist} {ifin} {cmd:,}
{opt prices(varlist)} {c |} {opt lnprices(varlist)}
{opt expenditure(varname)} {c |} {opt lnexpenditure(varname)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt :{opt prices(varlist)}}prices of goods in the system (same order as {it:sharelist}){p_end}
{synopt :{opt lnprices(varlist)}}log prices of goods in the system{p_end}
{synopt :{opt expenditure(varname)}}total system expenditure{p_end}
{synopt :{opt lnexpenditure(varname)}}log total system expenditure{p_end}
{synopt :{opt demographics(varlist)}}demographic variables for Ray (1983) scaling; required for censoring{p_end}
{synopt :{opt anot(#)}}constant in the translog price index; default is {cmd:anot(10.0)}{p_end}
{synopt :{opt noquadratic}}estimate linear AIDS instead of QUAIDS{p_end}
{synopt :{opt nocensor}}estimate uncensored demand system without Shonkwiler-Yen correction{p_end}

{syntab:Control function and selection design}
{synopt :{opt ivexp(varlist)}}excluded instruments for endogenous expenditure; internally builds the residual used in both model stages{p_end}
{synopt :{opt control_function(varname)}}externally generated reduced-form residual entering latent demand shares{p_end}
{synopt :{opt selection_control_function(varname)}}residual entering each first-stage Probit with equation-specific coefficients{p_end}
{synopt :{opt selection_prices(varlist)}}ordered subset of demand prices used by the first-stage Probits{p_end}
{synopt :{opt selection_noprices}}omit all price variables from the first-stage Probits{p_end}
{synopt :{opt selection_covariates(varlist)}}ordered first-stage covariates, independent of Ray demographics{p_end}
{synopt :{opt selection_nocovariates}}omit all covariates from the first-stage Probits{p_end}
{synopt :{opt selection_noexpenditure}}omit log expenditure from the first-stage Probits{p_end}

{syntab:Estimation & Optimizer}
{synopt :{opt method(method)}}estimation method: {cmd:ifgnls} (default), {cmd:fgnls}, or {cmd:nls}{p_end}
{synopt :{opt algorithm(alg)}}optimizer algorithm: {cmd:gn} (Gauss-Newton, default) or {cmd:lm} (Levenberg-Marquardt){p_end}
{synopt :{opt start(type)}}starting values: {cmd:zero} (default, matches Stata) or {cmd:linear} (linearized AIDS start){p_end}
{synopt :{opt initial(matname)}}row vector (matrix name) of initial free parameters for warm-starting; e.g. {cmd:e(b_est)} from a prior run{p_end}
{synopt :{opt sigma_initial(matname)}}initial residual covariance matrix (matrix name), used with {opt initial()} for warm-starting; e.g. {cmd:e(Sigma)} from a prior run{p_end}
{synopt :{opt vce_sigma(type)}}GLS covariance convention: {cmd:objective} (default) or {cmd:final}; NLS always uses its identity-weighted sandwich{p_end}
{synopt :{opt param_tol(#)}}inner relative parameter-change threshold; default is {cmd:param_tol(1e-5)}{p_end}
{synopt :{opt objective_tol(#)}}inner relative weighted-SSR-change threshold; default is {cmd:objective_tol(1e-7)}{p_end}
{synopt :{opt gn_tol(#)}}inner scaled Gauss-Newton threshold; default is {cmd:gn_tol(1e-5)}{p_end}
{synopt :{opt outer_param_tol(#)}}IFGNLS relative parameter-change threshold, required in two consecutive rounds; default is {cmd:outer_param_tol(1e-5)}{p_end}
{synopt :{opt max_iter(#)}}maximum inner Gauss-Newton iterations per stage; default is {cmd:max_iter(300)}{p_end}
{synopt :{opt max_outer(#)}}maximum numbered IFGNLS estimation stage, including initial NLS and FGNLS; default is {cmd:max_outer(200)}{p_end}
{synopt :{opt chunk(#)}}observation block size for accumulating normal equations; default is {cmd:chunk(2000)}{p_end}

{syntab:Bootstrap & Performance}
{synopt :{opt reps(#)}}number of bootstrap replications; default is {cmd:reps(0)} (disabled){p_end}
{synopt :{opt bootstrap_start(type)}}bootstrap starting values: {cmd:zero} (default) or {cmd:warm} (fast warm-start){p_end}
{synopt :{opt seed(#)}}random number seed for bootstrap{p_end}
{synopt :{opt n_jobs(#)}}number of parallel CPU cores for bootstrap; default is {cmd:n_jobs(1)}{p_end}
{synopt :{opt blas_threads(#)}}BLAS threads used by each estimation process; default is {cmd:blas_threads(1)}{p_end}
{synopt :{opt mp_context(method)}}Python multiprocessing start method; safe default is {cmd:spawn}{p_end}
{synopt :{opt rep_timeout(#)}}cooperative plus parent-watchdog time limit in seconds for each bootstrap replication; 0 disables it{p_end}
{synopt :{opt nolog}}suppress estimation iteration log{p_end}
{synopt :{opt gnlog}}print detailed step-by-step Gauss-Newton optimization logs{p_end}
{synopt :{opt level(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synoptline}
{pstd}
The inner solve stops when any one of its parameter-change, objective-change,
or scaled Gauss-Newton criteria passes. Parameter change uses the maximum
absolute step divided by 1 plus the absolute old parameter. Objective change
uses the weighted SSR decrease divided by the absolute old weighted SSR
(with denominator floor 1e-300). The GN criterion uses the absolute inner
product of the undamped GN direction and score, divided by that same SSR
scale. If no improving step is accepted, only the GN criterion can certify
convergence.{p_end}

{pstd}
The outer IFGNLS criterion is the maximum absolute parameter change divided
by the absolute old parameter plus 1e-8. It must be below {opt outer_param_tol()}
in two consecutive rounds, and the final inner solve must converge. The four
thresholds must be finite and positive. They stay constant during the fit and
are shared by every bootstrap replication, whether started from zero or warm
estimates. These are pyquaidsce's criteria, without a claim of exact equivalence
to Stata's native {cmd:nlsur} stopping rules.{p_end}

{pstd}
NLS parameter covariance uses A^-1 B A^-1, with A=sum J'J and
B=sum J'Sigma_hat J; Sigma_hat is the final residual covariance divided by N.
It assumes a common covariance across independent observations. FGNLS and
IFGNLS retain their existing GLS covariance convention.{p_end}



{marker description}{...}
{title:Description}

{pstd}
You can use either {cmd:pyquaidsce} or {cmd:quaidsce}. Both commands run the same estimator. If you already have the older Stata {cmd:quaidsce} installed and want to avoid a name conflict, use {cmd:pyquaidsce}.

{pstd}
{cmd:pyquaidsce} estimates the Quadratic Almost Ideal Demand System (QUAIDS) of Banks, Blundell, and Lewbel (1997) with Ray (1983) demographic scaling and the Shonkwiler & Yen (1999) two-step correction for zero budget shares.

{pstd}
It provides a fast Stata front end powered by the {cmd:pyquaidsce} Python computation engine, achieving up to a {bf:44.6x speedup} under IFGNLS in benchmark tests. Point estimation and optional bootstrap replications run in a background Python process while Stata polls for progress, so the Stata GUI remains responsive.

{pstd}
The Shonkwiler-Yen correction always uses the Probit linear index. With
{cmd:ivexp()}, log expenditure is regressed
on log prices, Ray demographics, the excluded instruments, and a constant. The
generated residual enters both the participation Probits and latent demand
equations, with equation-specific coefficients. The internal bootstrap
re-estimates this reduced form in every replication. Bootstrap remains disabled
for residuals supplied through {cmd:control_function()} or
{cmd:selection_control_function()}, because those residuals cannot be rebuilt
from the information supplied to the command.


{marker examples}{...}
{title:Examples}

{pstd}Load household consumption data and estimate a 4-good censored QUAIDS model with IFGNLS:{p_end}

{phang2}{cmd:. use mydata.dta, clear}{p_end}
{phang2}{cmd:. quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) method(ifgnls)}{p_end}

{pstd}Run estimation with 200 parallel bootstrap replications across 4 CPU cores:{p_end}

{phang2}{cmd:. quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) reps(200) n_jobs(4) mp_context(spawn) rep_timeout(900) seed(123456)}{p_end}

{pstd}Use an externally generated demand residual and a distinct selection design:{p_end}

{phang2}{cmd:. quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) control_function(vhat) selection_control_function(vhat_sel) selection_prices(p3 p1) selection_covariates(urban) selection_noexpenditure reps(0)}{p_end}

{pstd}Instrument endogenous expenditure internally and use a full bootstrap:{p_end}

{phang2}{cmd:. quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) ivexp(log_income employment) reps(200) seed(12345)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:pyquaidsce} stores the following in {cmd:e()}:

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(ll)}}log-likelihood{p_end}
{synopt:{cmd:e(anot)}}price index constant{p_end}
{synopt:{cmd:e(ndemo)}}number of demographic variables{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(reduced_form_r2)}}R-squared from the internal log-expenditure reduced form{p_end}
{synopt:{cmd:e(excluded_iv_F)}}classical joint F statistic for the excluded instruments{p_end}
{synopt:{cmd:e(excluded_iv_p)}}p-value for the excluded-instrument F test{p_end}

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:pyquaidsce}{p_end}
{synopt:{cmd:e(title)}}model title{p_end}
{synopt:{cmd:e(method)}}estimation method ({cmd:ifgnls}, {cmd:fgnls}, {cmd:nls}){p_end}
{synopt:{cmd:e(control_function)}}demand control-function variable, if supplied{p_end}
{synopt:{cmd:e(selection_control_function)}}selection control-function variable, if supplied{p_end}
{synopt:{cmd:e(ivexp)}}excluded expenditure instruments, if supplied{p_end}

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}full coefficient vector{p_end}
{synopt:{cmd:e(b_est)}}vector of free estimated structural parameters (can be passed to {cmd:initial()} for warm-starting in subsequent runs){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(Sigma)}}residual covariance matrix (can be passed to {cmd:sigma_initial()} for warm-starting in subsequent runs){p_end}
{synopt:{cmd:e(elas_i)}}expenditure (income) elasticities{p_end}
{synopt:{cmd:e(elas_u)}}uncompensated (Marshallian) price elasticities matrix{p_end}
{synopt:{cmd:e(elas_c)}}compensated (Hicksian) price elasticities matrix{p_end}
{synopt:{cmd:e(reduced_form_b)}}internal log-expenditure reduced-form coefficients{p_end}
{synopt:{cmd:e(reduced_form_V)}}classical reduced-form coefficient covariance matrix{p_end}


{marker author}{...}
{title:Author}

{pstd}
{bf:Sina Amiri}{break}
Department of Economics, Shiraz University, Shiraz, Iran{break}
Email: {browse "mailto:sinaamiri9000@gmail.com":sinaamiri9000@gmail.com}{break}
GitHub: {browse "https://github.com/sinaamiri9000-collab/pyquaidsce"}
