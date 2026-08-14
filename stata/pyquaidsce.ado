*! version 1.0.1  14aug2026
*! pyquaidsce: Censored QUAIDS demand system estimation in Stata using Python engine
*! Author: Sina Amiri (Department of Economics, Shiraz University)

program define pyquaidsce, eclass
    version 16.0

    syntax varlist(min=3 numeric) [if] [in], ///
        [ prices(varlist numeric) lnprices(varlist numeric) ///
          expenditure(varname numeric) lnexpenditure(varname numeric) ///
          demographics(varlist numeric) anot(real 10.0) ///
          method(string) algorithm(string) reps(integer 0) ///
          first_stage_predict(string) strict_stata(string) ///
          seed(integer -1) n_jobs(integer 1) noquadratic nocensor nolog level(cilevel) ]

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

    // 2. Mark estimation sample
    marksample touse
    markout `touse' `p_vars' `exp_var' `demographics'

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

    if "`algorithm'" == "" local algorithm "gn"
    local algorithm = lower("`algorithm'")
    if !inlist("`algorithm'", "gn", "lm") {
        display as error "algorithm must be 'gn' or 'lm'"
        exit 198
    }

    if "`first_stage_predict'" == "" local first_stage_predict "pr"
    local first_stage_predict = lower("`first_stage_predict'")
    if !inlist("`first_stage_predict'", "pr", "xb") {
        display as error "first_stage_predict must be 'pr' or 'xb'"
        exit 198
    }

    local is_quad = ("`quadratic'" == "")
    local is_censor = ("`censor'" == "")
    local is_strict = ("`strict_stata'" != "false" & "`strict_stata'" != "0")
    local is_verbose = ("`log'" == "")

    // 4. Temporary matrices for ereturn
    tempname b V elas_i elas_u elas_c

    // 5. Call Python Bridge
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

    python: from pyquaidsce.stata_bridge import run_from_stata; import sfi; run_from_stata(shares_str=sfi.Macro.getLocal("varlist"), prices_str=sfi.Macro.getLocal("p_vars"), expenditure_str=sfi.Macro.getLocal("exp_var"), demographics_str=sfi.Macro.getLocal("demographics"), anot=float(sfi.Macro.getLocal("anot")), method=sfi.Macro.getLocal("method"), algorithm=sfi.Macro.getLocal("algorithm"), reps=int(sfi.Macro.getLocal("reps")), seed=int(sfi.Macro.getLocal("seed")), n_jobs=int(sfi.Macro.getLocal("n_jobs")), first_stage_predict=sfi.Macro.getLocal("first_stage_predict"), strict_stata=bool(int(sfi.Macro.getLocal("is_strict"))), quadratic=bool(int(sfi.Macro.getLocal("is_quad"))), censor=bool(int(sfi.Macro.getLocal("is_censor"))), is_lnprices=bool(int(sfi.Macro.getLocal("is_lnp"))), is_lnexp=bool(int(sfi.Macro.getLocal("is_lnexp"))), verbose=bool(int(sfi.Macro.getLocal("is_verbose"))), b_mat_name=sfi.Macro.getLocal("b"), v_mat_name=sfi.Macro.getLocal("V"), elas_i_name=sfi.Macro.getLocal("elas_i"), elas_u_name=sfi.Macro.getLocal("elas_u"), elas_c_name=sfi.Macro.getLocal("elas_c"), touse_var=sfi.Macro.getLocal("touse"))

    // 6. Post estimation results to e()
    ereturn post `b' `V', esample(`touse')
    ereturn scalar N = scalar(r_nobs)
    ereturn scalar ll = scalar(r_llf)
    ereturn scalar anot = scalar(r_anot)
    ereturn scalar ndemo = scalar(r_ndemo)
    ereturn scalar converged = scalar(r_converged)
    ereturn scalar n_outer = scalar(r_n_outer)
    ereturn scalar n_gn = scalar(r_n_gn)

    ereturn matrix elas_i = `elas_i'
    ereturn matrix elas_u = `elas_u'
    ereturn matrix elas_c = `elas_c'

    ereturn local cmd "pyquaidsce"
    ereturn local cmdline "pyquaidsce `0'"
    ereturn local title "`model_title'"
    ereturn local method "`method'"
    ereturn local predict "`first_stage_predict'"
    ereturn local shares "`varlist'"
    ereturn local prices "`p_vars'"
    ereturn local demographics "`demographics'"

    // 7. Display Stata output
    display _n as text "`model_title'"
    display as text "{hline 78}"
    display as text "Number of obs          =" as result %10.0f e(N)
    display as text "Number of demographics =" as result %10.0f e(ndemo)
    display as text "Alpha_0                =" as result %10.4f e(anot)
    display as text "Log-likelihood         =" as result %10.4f e(ll)
    if `reps' > 0 {
        display as text "Bootstrap replications =" as result %10.0f `reps'
    }

    ereturn display, level(`level')

    display _n as text "Expenditure (income) elasticities, at means:"
    matrix list e(elas_i), noheader format(%10.6f)

    display _n as text "Uncompensated (Marshallian) price elasticities [row = good, column = price]:"
    matrix list e(elas_u), noheader format(%10.6f)

    display _n as text "Compensated (Hicksian) price elasticities [row = good, column = price]:"
    matrix list e(elas_c), noheader format(%10.6f)

end
