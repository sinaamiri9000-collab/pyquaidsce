{smcl}
{* *! version 1.1.0  14aug2026}{...}
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

{syntab:Estimation & Optimizer}
{synopt :{opt method(method)}}estimation method: {cmd:ifgnls} (default), {cmd:fgnls}, or {cmd:nls}{p_end}
{synopt :{opt algorithm(alg)}}optimizer algorithm: {cmd:gn} (Gauss-Newton, default) or {cmd:lm} (Levenberg-Marquardt){p_end}
{synopt :{opt first_stage_predict(type)}}{cmd:pr} (default, matches Stata) or {cmd:xb} (textbook linear index){p_end}
{synopt :{opt strict_stata(bool)}}{cmd:true} (default, matches Stata formulas) or {cmd:false} (corrected formulas){p_end}

{syntab:Bootstrap & Performance}
{synopt :{opt reps(#)}}number of bootstrap replications; default is {cmd:reps(0)} (disabled){p_end}
{synopt :{opt seed(#)}}random number seed for bootstrap{p_end}
{synopt :{opt n_jobs(#)}}number of parallel CPU cores for bootstrap; default is {cmd:n_jobs(1)}{p_end}
{synopt :{opt nolog}}suppress estimation iteration log{p_end}
{synopt :{opt level(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:pyquaidsce} estimates the Quadratic Almost Ideal Demand System (QUAIDS) of Banks, Blundell, and Lewbel (1997) with Ray (1983) demographic scaling and the Shonkwiler & Yen (1999) two-step correction for zero budget shares.

{pstd}
It provides a fast Stata front end powered by the {cmd:pyquaidsce} Python computation engine, achieving up to a {bf:44.6x speedup} under IFGNLS in benchmark tests.


{marker examples}{...}
{title:Examples}

{pstd}Load household consumption data and estimate a 4-good censored QUAIDS model with IFGNLS:{p_end}

{phang2}{cmd:. use mydata.dta, clear}{p_end}
{phang2}{cmd:. pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) method(ifgnls)}{p_end}

{pstd}Run estimation with 200 parallel bootstrap replications across 4 CPU cores:{p_end}

{phang2}{cmd:. pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) reps(200) n_jobs(4) seed(123456)}{p_end}


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

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(elas_i)}}expenditure (income) elasticities{p_end}
{synopt:{cmd:e(elas_u)}}uncompensated (Marshallian) price elasticities matrix{p_end}
{synopt:{cmd:e(elas_c)}}compensated (Hicksian) price elasticities matrix{p_end}


{marker author}{...}
{title:Author}

{pstd}
{bf:Sina Amiri}{break}
Department of Economics, Shiraz University, Shiraz, Iran{break}
GitHub: {browse "https://github.com/sinaamiri9000-collab/pyquaidsce"}
