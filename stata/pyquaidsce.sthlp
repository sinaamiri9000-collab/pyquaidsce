{smcl}
{* *! version 1.3.0  22aug2026}{...}
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
{synopt :{opt control_function(varname)}}externally generated reduced-form residual entering latent demand shares{p_end}
{synopt :{opt selection_control_function(varname)}}residual entering each first-stage Probit with equation-specific coefficients{p_end}
{synopt :{opt selection_prices(varlist)}}ordered subset of demand prices used by the first-stage Probits{p_end}
{synopt :{opt selection_noprices}}omit all price variables from the first-stage Probits{p_end}
{synopt :{opt selection_covariates(varlist)}}ordered first-stage covariates, independent of Ray demographics{p_end}
{synopt :{opt selection_nocovariates}}omit all covariates from the first-stage Probits{p_end}
{synopt :{opt selection_noexpenditure}}omit log expenditure from the first-stage Probits{p_end}

{syntab:Estimation & Optimizer}
{synopt :{opt method(method)}}estimation method: {cmd:ifgnls} (default), {cmd:fgnls}, or {cmd:nls}{p_end}
{synopt :{opt stop_rule(rule)}}stopping rule: {cmd:standard} (default, disjunctive rule matching Stata's {opt tolerance}/{opt ltolerance}/{opt nrtolerance}) or {cmd:tight} (strict scaled-gradient rule){p_end}
{synopt :{opt algorithm(alg)}}optimizer algorithm: {cmd:gn} (Gauss-Newton, default) or {cmd:lm} (Levenberg-Marquardt){p_end}
{synopt :{opt start(type)}}starting values: {cmd:zero} (default, matches Stata) or {cmd:linear} (linearized AIDS start){p_end}
{synopt :{opt initial(matname)}}row vector (matrix name) of initial free parameters for warm-starting; e.g. {cmd:e(b_est)} from a prior run{p_end}
{synopt :{opt sigma_initial(matname)}}initial residual covariance matrix (matrix name), used with {opt initial()} for warm-starting; e.g. {cmd:e(Sigma)} from a prior run{p_end}
{synopt :{opt vce_sigma(type)}}covariance used in the S.E. formula: {cmd:objective} (default, matches Stata) or {cmd:final}{p_end}
{synopt :{opt tol(#)}}objective relative-change tolerance; default is {cmd:tol(1e-13)}{p_end}
{synopt :{opt nrtol_stop(#)}}scaled relative gradient stopping tolerance; default is {cmd:nrtol_stop(1e-12)}{p_end}
{synopt :{opt sigma_tol(#)}}outer fixed-point parameter tolerance for IFGNLS; default is {cmd:sigma_tol(1e-11)}{p_end}
{synopt :{opt inner_nrtol_early(#)}}early-stage inner Gauss-Newton tolerance during inexact-outer IFGNLS; default is {cmd:inner_nrtol_early(1e-8)}{p_end}
{synopt :{opt max_iter(#)}}maximum inner Gauss-Newton iterations per stage; default is {cmd:max_iter(300)}{p_end}
{synopt :{opt max_outer(#)}}maximum outer covariance updates for IFGNLS; default is {cmd:max_outer(200)}{p_end}
{synopt :{opt chunk(#)}}observation block size for accumulating normal equations; default is {cmd:chunk(2000)}{p_end}
{synopt :{opt first_stage_predict(type)}}{cmd:xb} (default, textbook linear index) or {cmd:pr} (legacy Stata-style prediction){p_end}
{synopt :{opt strict_stata(bool)}}{cmd:false} (default, corrected textbook formulas) or {cmd:true} (reproduces the original ado's elasticity formulas){p_end}

{syntab:Bootstrap & Performance}
{synopt :{opt reps(#)}}number of bootstrap replications; default is {cmd:reps(0)} (disabled){p_end}
{synopt :{opt bootstrap_start(type)}}bootstrap starting values: {cmd:zero} (default) or {cmd:warm} (fast warm-start){p_end}
{synopt :{opt boot_sigma_tol(#)}}outer covariance tolerance inside each bootstrap replication; default is {cmd:boot_sigma_tol(1e-7)}{p_end}
{synopt :{opt seed(#)}}random number seed for bootstrap{p_end}
{synopt :{opt n_jobs(#)}}number of parallel CPU cores for bootstrap; default is {cmd:n_jobs(1)}{p_end}
{synopt :{opt mp_context(method)}}Python multiprocessing start method; safe default is {cmd:spawn}{p_end}
{synopt :{opt rep_timeout(#)}}cooperative plus parent-watchdog time limit in seconds for each bootstrap replication; 0 disables it{p_end}
{synopt :{opt nolog}}suppress estimation iteration log{p_end}
{synopt :{opt gnlog}}print detailed step-by-step Gauss-Newton optimization logs{p_end}
{synopt :{opt level(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:pyquaidsce} estimates the Quadratic Almost Ideal Demand System (QUAIDS) of Banks, Blundell, and Lewbel (1997) with Ray (1983) demographic scaling and the Shonkwiler & Yen (1999) two-step correction for zero budget shares.

{pstd}
It provides a fast Stata front end powered by the {cmd:pyquaidsce} Python computation engine, achieving up to a {bf:44.6x speedup} under IFGNLS in benchmark tests. Point estimation and optional bootstrap replications run in a background Python process while Stata polls for progress, so the Stata GUI remains responsive.

{pstd}
Control-function and custom-selection options require
{cmd:first_stage_predict(xb)} and are currently restricted to
{cmd:reps(0)}. A valid bootstrap for a generated residual must re-estimate its
reduced form inside every replication; passing a precomputed residual unchanged
would understate uncertainty.


{marker examples}{...}
{title:Examples}

{pstd}Load household consumption data and estimate a 4-good censored QUAIDS model with IFGNLS:{p_end}

{phang2}{cmd:. use mydata.dta, clear}{p_end}
{phang2}{cmd:. pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) method(ifgnls)}{p_end}

{pstd}Run estimation with 200 parallel bootstrap replications across 4 CPU cores:{p_end}

{phang2}{cmd:. pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) reps(200) n_jobs(4) mp_context(spawn) rep_timeout(900) seed(123456)}{p_end}

{pstd}Use an externally generated demand residual and a distinct selection design:{p_end}

{phang2}{cmd:. pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) control_function(vhat) selection_control_function(vhat_sel) selection_prices(p3 p1) selection_covariates(urban) selection_noexpenditure first_stage_predict(xb) reps(0)}{p_end}


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

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:pyquaidsce}{p_end}
{synopt:{cmd:e(title)}}model title{p_end}
{synopt:{cmd:e(method)}}estimation method ({cmd:ifgnls}, {cmd:fgnls}, {cmd:nls}){p_end}
{synopt:{cmd:e(control_function)}}demand control-function variable, if supplied{p_end}
{synopt:{cmd:e(selection_control_function)}}selection control-function variable, if supplied{p_end}

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}full coefficient vector{p_end}
{synopt:{cmd:e(b_est)}}vector of free estimated structural parameters (can be passed to {cmd:initial()} for warm-starting in subsequent runs){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(Sigma)}}residual covariance matrix (can be passed to {cmd:sigma_initial()} for warm-starting in subsequent runs){p_end}
{synopt:{cmd:e(elas_i)}}expenditure (income) elasticities{p_end}
{synopt:{cmd:e(elas_u)}}uncompensated (Marshallian) price elasticities matrix{p_end}
{synopt:{cmd:e(elas_c)}}compensated (Hicksian) price elasticities matrix{p_end}


{marker author}{...}
{title:Author}

{pstd}
{bf:Sina Amiri}{break}
Department of Economics, Shiraz University, Shiraz, Iran{break}
Email: {browse "mailto:sinaamiri9000@gmail.com":sinaamiri9000@gmail.com}{break}
GitHub: {browse "https://github.com/sinaamiri9000-collab/pyquaidsce"}
