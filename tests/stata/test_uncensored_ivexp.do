* Run from the repository root with the matching Python backend installed.
version 16.0
clear all
set more off
adopath ++ "stata"
use "tests/data/small4.dta", clear

* Artificial instruments exercise the interface, not economic identification.
set seed 918
generate double test_iv = ln(total) + rnormal(0, 0.5)
generate double test_iv_sq = test_iv^2

foreach specification in aids quaids {
    local quadratic_option ""
    if "`specification'" == "aids" local quadratic_option "noquadratic"
    pyquaidsce sw1 sw2 sw4 sw9, prices(p1 p2 p4 p9) expenditure(total) ///
        demographics(x1 x2) ivexp(test_iv test_iv_sq) nocensor ///
        `quadratic_option' anot(10) reps(0) nolog
    assert e(N) == 2000
    assert e(converged) == 1
    assert abs(_b[cfcoef:cfcoef_1] + _b[cfcoef:cfcoef_2] + ///
        _b[cfcoef:cfcoef_3] + _b[cfcoef:cfcoef_4]) < 1e-12
    assert e(excluded_iv_F) >= 0
    tempname EU EC
    matrix `EU' = e(elas_u)
    matrix `EC' = e(elas_c)
    forvalues i = 1/4 {
        forvalues j = 1/4 {
            assert abs(_b[ELAS_UNCOMP:e_`i'_`j'] - el(`EU', `i', `j')) < 1e-12
            assert abs(_b[ELAS_COMP:e_`i'_`j'] - el(`EC', `i', `j')) < 1e-12
        }
    }
}

* The common command name must forward the same options, including bootstrap.
quaidsce sw1 sw2 sw4 sw9, prices(p1 p2 p4 p9) expenditure(total) ///
    demographics(x1 x2) ivexp(test_iv test_iv_sq) nocensor ///
    anot(10) reps(4) seed(573) n_jobs(1) bootstrap_start(warm) nolog
assert e(converged) == 1
assert e(boot_reps) == 4
assert _se[cfcoef:cfcoef_4] > 0
assert _se[ELAS_INC:e_1] > 0

* Participation options still require censoring correction.
capture pyquaidsce sw1 sw2 sw4 sw9, prices(p1 p2 p4 p9) ///
    expenditure(total) demographics(x1 x2) nocensor selection_noexpenditure
assert _rc == 198

* A precomputed residual still cannot use the internal bootstrap.
generate double residual = test_iv - ln(total)
capture pyquaidsce sw1 sw2 sw4 sw9, prices(p1 p2 p4 p9) ///
    expenditure(total) demographics(x1 x2) nocensor control_function(residual) reps(4)
assert _rc == 198
