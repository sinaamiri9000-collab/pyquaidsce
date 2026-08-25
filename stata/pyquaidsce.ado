*! version 1.5.0  25aug2026
*! pyquaidsce: Censored QUAIDS demand system estimation in Stata using Python engine
*! Author: Sina Amiri (Department of Economics, Shiraz University)

program define pyquaidsce, eclass
    version 16.0

    syntax varlist(min=3 numeric) [if] [in], ///
        [ prices(varlist numeric) lnprices(varlist numeric) ///
          expenditure(varname numeric) lnexpenditure(varname numeric) ///
          demographics(varlist numeric) ivexp(varlist numeric) anot(real 10.0) ///
          control_function(varname numeric) ///
          selection_control_function(varname numeric) ///
          selection_prices(varlist numeric) selection_noprices ///
          selection_covariates(varlist numeric) selection_nocovariates ///
          selection_noexpenditure ///
          method(string) algorithm(string) start(string) reps(integer 0) ///
          stop_rule(string) bootstrap_start(string) ///
          first_stage_predict(string) strict_stata(string) ///
          vce_sigma(string) ///
          initial(string) sigma_initial(string) ///
          tol(real 1e-13) max_outer(integer 200) max_iter(integer 300) ///
          chunk(integer 2000) nrtol_stop(real 1e-12) ///
          inner_nrtol_early(real 1e-8) sigma_tol(real 1e-11) ///
          boot_sigma_tol(real 1e-7) ///
          seed(integer -1) n_jobs(integer 1) mp_context(string) ///
          rep_timeout(real 0) noquadratic nocensor nolog gnlog level(cilevel) ]

    // 1. Validate inputs
    if "`prices'" == "" & "`lnprices'" == "" {
        display as error "must specify either prices() or lnprices()"
        exit 198
    }
    if "`prices'" != "" & "`lnprices'" != "" {
        display as error "cannot specify both prices() and lnprices()"
        exit 198
    }
    if "`expenditure'" == "" & "`lnexpenditure'" == "" {
        display as error "must specify either expenditure() or lnexpenditure()"
        exit 198
    }
    if "`expenditure'" != "" & "`lnexpenditure'" != "" {
        display as error "cannot specify both expenditure() and lnexpenditure()"
        exit 198
    }
    if "`ivexp'" != "" & ///
       ("`control_function'" != "" | "`selection_control_function'" != "") {
        display as error "ivexp() cannot be combined with control_function() or selection_control_function()"
        exit 198
    }

    local p_vars = cond("`prices'" != "", "`prices'", "`lnprices'")
    local exp_var = cond("`expenditure'" != "", "`expenditure'", "`lnexpenditure'")
    local is_lnp = ("`lnprices'" != "")
    local is_lnexp = ("`lnexpenditure'" != "")

    local n_shares : word count `varlist'
    local n_prices : word count `p_vars'
    if `n_shares' != `n_prices' {
        display as error "number of prices (`n_prices') must equal number of budget shares (`n_shares')"
        exit 198
    }
    if "`selection_prices'" != "" & "`selection_noprices'" != "" {
        display as error "cannot combine selection_prices() and selection_noprices"
        exit 198
    }
    if "`selection_covariates'" != "" & "`selection_nocovariates'" != "" {
        display as error "cannot combine selection_covariates() and selection_nocovariates"
        exit 198
    }

    // 2. Mark estimation sample
    marksample touse
    markout `touse' `p_vars' `exp_var' `demographics' `ivexp' ///
        `control_function' `selection_control_function' ///
        `selection_prices' `selection_covariates'

    quietly count if `touse'
    if r(N) == 0 {
        error 2000
    }
    local N_obs = r(N)

    // 3. Set default options
    if "`method'" == "" local method "ifgnls"
    local method = lower("`method'")
    if !inlist("`method'", "nls", "fgnls", "ifgnls") {
        display as error "method must be 'nls', 'fgnls', or 'ifgnls'"
        exit 198
    }

    if "`stop_rule'" == "" local stop_rule "standard"
    local stop_rule = lower("`stop_rule'")
    if !inlist("`stop_rule'", "standard", "tight") {
        display as error "stop_rule must be 'standard' or 'tight'"
        exit 198
    }

    if "`bootstrap_start'" == "" local bootstrap_start "zero"
    local bootstrap_start = lower("`bootstrap_start'")
    if !inlist("`bootstrap_start'", "zero", "warm") {
        display as error "bootstrap_start must be 'zero' or 'warm'"
        exit 198
    }

    if "`algorithm'" == "" local algorithm "gn"
    local algorithm = lower("`algorithm'")
    if !inlist("`algorithm'", "gn", "lm") {
        display as error "algorithm must be 'gn' or 'lm'"
        exit 198
    }

    if "`start'" == "" local start "zero"
    local start = lower("`start'")
    if !inlist("`start'", "zero", "linear") {
        display as error "start must be 'zero' or 'linear'"
        exit 198
    }

    if "`first_stage_predict'" == "" local first_stage_predict "xb"
    local first_stage_predict = lower("`first_stage_predict'")
    if !inlist("`first_stage_predict'", "pr", "xb") {
        display as error "first_stage_predict must be 'pr' or 'xb'"
        exit 198
    }

    if "`vce_sigma'" == "" local vce_sigma "objective"
    local vce_sigma = lower("`vce_sigma'")
    if !inlist("`vce_sigma'", "objective", "final") {
        display as error "vce_sigma must be 'objective' or 'final'"
        exit 198
    }

    // initial()/sigma_initial() accept a Stata matrix name (e.g. from a prior
    // run's e(b)/e(V) or any k x 1 vector / m x m matrix in memory).
    if "`initial'" != "" {
        capture confirm matrix `initial'
        if _rc {
            display as error "initial() must be the name of an existing Stata matrix"
            exit 198
        }
    }
    if "`sigma_initial'" != "" {
        capture confirm matrix `sigma_initial'
        if _rc {
            display as error "sigma_initial() must be the name of an existing Stata matrix"
            exit 198
        }
    }
    if `tol' <= 0 | `nrtol_stop' <= 0 | `sigma_tol' <= 0 ///
        | `boot_sigma_tol' <= 0 {
        display as error "convergence tolerances must be positive"
        exit 198
    }
    if `max_outer' < 2 | `max_iter' < 1 | `chunk' < 1 {
        display as error "max_outer must be >= 2; max_iter and chunk must be >= 1"
        exit 198
    }
    if `inner_nrtol_early' <= 0 {
        display as error "inner_nrtol_early must be positive"
        exit 198
    }

    local is_quad = ("`quadratic'" == "")
    local is_censor = ("`censor'" == "")
    // Default is now strict_stata(false): corrected/textbook formulas.
    if "`strict_stata'" != "" {
        local strict_stata = lower("`strict_stata'")
        if !inlist("`strict_stata'", "true", "false", "1", "0") {
            display as error "strict_stata must be 'true' or 'false'"
            exit 198
        }
    }
    local is_strict = ("`strict_stata'" == "true" | "`strict_stata'" == "1")
    local is_verbose = ("`log'" == "")
    local is_gn_verbose = ("`gnlog'" != "")
    local selection_prices_specified = ///
        ("`selection_prices'" != "" | "`selection_noprices'" != "")
    local selection_covariates_specified = ///
        ("`selection_covariates'" != "" | "`selection_nocovariates'" != "")
    local selection_expenditure_on = ("`selection_noexpenditure'" == "")
    local external_cf_active = ///
        ("`control_function'" != "" | "`selection_control_function'" != "")
    local extension_active = ///
        (`external_cf_active' | "`ivexp'" != "" | ///
         `selection_prices_specified' | `selection_covariates_specified' | ///
         !`selection_expenditure_on')
    if `extension_active' & !`is_censor' {
        display as error "control-function/selection extensions require censoring"
        exit 198
    }
    if `extension_active' & "`first_stage_predict'" != "xb" {
        display as error "control-function/selection extensions require first_stage_predict(xb)"
        exit 198
    }
    if `external_cf_active' & `reps' > 0 {
        display as error "bootstrap is disabled for precomputed control residuals"
        display as text  "Use ivexp() so the reduced form is rebuilt inside each bootstrap replication."
        exit 198
    }
    if `rep_timeout' < 0 {
        display as error "rep_timeout must be nonnegative (0 disables it)"
        exit 198
    }

    // 4. Temporary matrices for ereturn
    tempname b V elas_i elas_u elas_c b_est Sigma rf_b rf_V

    // 5. Check Python package
    capture python: import pyquaidsce
    if _rc != 0 {
        display as error "The 'pyquaidsce' Python package is not found in Stata's Python environment."
        display as text  "Attempting automatic installation via pip..."
        capture python: import subprocess, sys; subprocess.check_call([sys.executable, "-m", "pip", "install", "pyquaidsce"])
        if _rc != 0 {
            display as error "Automatic installation failed. Please install pyquaidsce manually by running:"
            display as text  "   pip install pyquaidsce"
            exit 198
        }
    }

    // 6. Launch the complete estimation outside Stata's GUI process
    python: from pyquaidsce.stata_bridge import launch_from_stata, poll_bootstrap, load_stata_results, kill_bootstrap; import sfi; _rt=float(sfi.Macro.getLocal("rep_timeout")); launch_from_stata(shares_str=sfi.Macro.getLocal("varlist"), prices_str=sfi.Macro.getLocal("p_vars"), expenditure_str=sfi.Macro.getLocal("exp_var"), demographics_str=sfi.Macro.getLocal("demographics"), anot=float(sfi.Macro.getLocal("anot")), method=sfi.Macro.getLocal("method"), algorithm=sfi.Macro.getLocal("algorithm"), start=sfi.Macro.getLocal("start"), reps=int(sfi.Macro.getLocal("reps")), stop_rule=sfi.Macro.getLocal("stop_rule"), bootstrap_start=sfi.Macro.getLocal("bootstrap_start"), seed=int(sfi.Macro.getLocal("seed")), n_jobs=int(sfi.Macro.getLocal("n_jobs")), mp_context=sfi.Macro.getLocal("mp_context") or None, rep_timeout=_rt if _rt > 0 else None, first_stage_predict=sfi.Macro.getLocal("first_stage_predict"), strict_stata=bool(int(sfi.Macro.getLocal("is_strict"))), quadratic=bool(int(sfi.Macro.getLocal("is_quad"))), censor=bool(int(sfi.Macro.getLocal("is_censor"))), is_lnprices=bool(int(sfi.Macro.getLocal("is_lnp"))), is_lnexp=bool(int(sfi.Macro.getLocal("is_lnexp"))), ivexp_str=sfi.Macro.getLocal("ivexp"), control_function=sfi.Macro.getLocal("control_function"), selection_control_function=sfi.Macro.getLocal("selection_control_function"), selection_prices_str=sfi.Macro.getLocal("selection_prices"), selection_prices_specified=bool(int(sfi.Macro.getLocal("selection_prices_specified"))), selection_covariates_str=sfi.Macro.getLocal("selection_covariates"), selection_covariates_specified=bool(int(sfi.Macro.getLocal("selection_covariates_specified"))), selection_expenditure=bool(int(sfi.Macro.getLocal("selection_expenditure_on"))), verbose=bool(int(sfi.Macro.getLocal("is_verbose"))), vce_sigma=sfi.Macro.getLocal("vce_sigma"), initial_mat_name=sfi.Macro.getLocal("initial"), sigma_initial_mat_name=sfi.Macro.getLocal("sigma_initial"), tol=float(sfi.Macro.getLocal("tol")), max_outer=int(float(sfi.Macro.getLocal("max_outer"))), max_iter=int(float(sfi.Macro.getLocal("max_iter"))), chunk=int(float(sfi.Macro.getLocal("chunk"))), nrtol_stop=float(sfi.Macro.getLocal("nrtol_stop")), inner_nrtol_early=float(sfi.Macro.getLocal("inner_nrtol_early")), sigma_tol=float(sfi.Macro.getLocal("sigma_tol")), boot_sigma_tol=float(sfi.Macro.getLocal("boot_sigma_tol")), gn_verbose=bool(int(sfi.Macro.getLocal("is_gn_verbose"))), touse_var=sfi.Macro.getLocal("touse"))

    local _pyq_boot_done = 0
    local _pyq_boot_err = 0
    local _pyq_last_msg ""
    while !`_pyq_boot_done' & !`_pyq_boot_err' {
        sleep 500
        local _pyq_boot_msg ""
        python: poll_bootstrap()
        local _pyq_boot_done = scalar(_pyq_boot_done)
        local _pyq_boot_err = scalar(_pyq_boot_err)
        if "`_pyq_boot_msg'" != "" & "`_pyq_boot_msg'" != "`_pyq_last_msg'" {
            display as text "  `_pyq_boot_msg'"
            local _pyq_last_msg "`_pyq_boot_msg'"
        }
    }

    if `_pyq_boot_err' {
        display as error "Estimation failed: `_pyq_boot_errmsg'"
        capture python: kill_bootstrap()
        exit 498
    }

    python: load_stata_results(sfi.Macro.getLocal("b"), sfi.Macro.getLocal("V"), sfi.Macro.getLocal("elas_i"), sfi.Macro.getLocal("elas_u"), sfi.Macro.getLocal("elas_c"), sfi.Macro.getLocal("b_est"), sfi.Macro.getLocal("Sigma"), sfi.Macro.getLocal("rf_b"), sfi.Macro.getLocal("rf_V"))

    // 7. Post completed point-estimate and optional bootstrap results
    ereturn post `b' `V', esample(`touse')
    ereturn scalar N = scalar(r_nobs)
    ereturn scalar ll = scalar(r_llf)
    ereturn scalar anot = scalar(r_anot)
    ereturn scalar ndemo = scalar(r_ndemo)
    ereturn scalar converged = scalar(r_converged)
    ereturn scalar n_outer = scalar(r_n_outer)
    ereturn scalar n_gn = scalar(r_n_gn)
    if `reps' > 0 {
        ereturn scalar boot_reps = scalar(r_boot_reps_ok)
    }

    ereturn matrix elas_i = `elas_i'
    ereturn matrix elas_u = `elas_u'
    ereturn matrix elas_c = `elas_c'
    ereturn matrix b_est = `b_est'
    ereturn matrix Sigma = `Sigma'
    if "`ivexp'" != "" {
        ereturn matrix reduced_form_b = `rf_b'
        ereturn matrix reduced_form_V = `rf_V'
        ereturn scalar reduced_form_r2 = scalar(r_rf_r2)
        ereturn scalar excluded_iv_F = scalar(r_rf_F)
        ereturn scalar excluded_iv_p = scalar(r_rf_p)
        ereturn scalar excluded_iv_df1 = scalar(r_rf_df1)
        ereturn scalar excluded_iv_df2 = scalar(r_rf_df2)
    }

    ereturn local cmd "pyquaidsce"
    ereturn local cmdline "pyquaidsce `0'"
    ereturn local title "`model_title'"
    ereturn local method "`method'"
    ereturn local predict "`first_stage_predict'"
    ereturn local shares "`varlist'"
    ereturn local prices "`p_vars'"
    ereturn local demographics "`demographics'"
    ereturn local ivexp "`ivexp'"
    ereturn local control_function "`control_function'"
    ereturn local selection_control_function "`selection_control_function'"
    ereturn local selection_prices "`selection_prices'"
    ereturn local selection_covariates "`selection_covariates'"

    // 7. Display Stata output
    display _n as text "`model_title'"
    display as text "{hline 78}"
    display as text "Number of obs          =" as result %10.0f e(N)
    display as text "Number of demographics =" as result %10.0f e(ndemo)
    display as text "Alpha_0                =" as result %10.4f e(anot)
    display as text "Log-likelihood         =" as result %10.4f e(ll)
    if `reps' > 0 {
        display as text "Bootstrap replications =" as result %10.0f e(boot_reps) as text "/" as result %5.0f `reps'
    }
    if "`ivexp'" != "" {
        display as text "Reduced-form R-squared =" as result %10.4f e(reduced_form_r2)
        display as text "Excluded-instrument F =" as result %10.4f e(excluded_iv_F) ///
            as text "  F(" as result %3.0f e(excluded_iv_df1) as text "," ///
            as result %6.0f e(excluded_iv_df2) as text "), p=" ///
            as result %9.5f e(excluded_iv_p)
    }

    ereturn display, level(`level')

    display _n as text "Expenditure (income) elasticities, at means:"
    matrix list e(elas_i), noheader format(%10.6f)

    display _n as text "Uncompensated (Marshallian) price elasticities [row = good, column = price]:"
    matrix list e(elas_u), noheader format(%10.6f)

    display _n as text "Compensated (Hicksian) price elasticities [row = good, column = price]:"
    matrix list e(elas_c), noheader format(%10.6f)

end
