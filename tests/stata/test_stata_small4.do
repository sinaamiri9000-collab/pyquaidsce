* Simple Stata smoke test for the current checkout.
* First install the Python package OUTSIDE Stata, and install the Stata ado.
* Restart Stata before executing this file.
* No sys.executable overrides, pip calls, or Python path changes are made here.

version 16.0
clear all
set more off
adopath ++ "stata"

* Set to 1 only AFTER the point-estimate test passes.
local test_bootstrap 0

capture log close pub_smoke
log using "installed_small4_smoke.log", replace text name(pub_smoke)

display as text _newline "=== CHECK PYQUAIDSCE INSTALLATION ==="
capture noisily which pyquaidsce
if _rc {
    display as error "pyquaidsce.ado is not installed. Install it from the current checkout first."
    log close pub_smoke
    exit 111
}

* Check that the subprocess launcher is available in the installed backend.
* This import is quick: it does not install anything or estimate a model.
capture noisily python: from pyquaidsce.stata_bridge import _get_stata_python_executable; print("PASS: Python launcher is available")
if _rc {
    display as error "Matching Python backend not found. Install the current checkout in the Python selected by Stata, then restart Stata."
    log close pub_smoke
    exit 198
}

display as text _newline "=== LOAD LOCAL small4.dta ==="
use "tests/data/small4.dta", clear
assert _N == 2000

display as text _newline "=== IFGNLS POINT ESTIMATION: reps(0) ==="
timer clear 1
timer on 1
capture noisily pyquaidsce sw1 sw2 sw4 sw9, ///
    prices(p1 p2 p4 p9) expenditure(total) ///
    demographics(x1 x2) anot(10) method(ifgnls) ///
    reps(0) nolog
local rc = _rc
timer off 1
if `rc' != 0 {
    display as error "FAIL: pyquaidsce returned Stata code `rc'."
    log close pub_smoke
    exit `rc'
}
assert e(N) == 2000
assert e(converged) == 1
scalar ll_point = e(ll)
display as result "PASS: Stata launched the installed Python worker and received results."
display as result "N = " e(N) "  Log-likelihood = " %16.6f e(ll)
timer list 1

tempname EU EC
matrix `EU' = e(elas_u)
matrix `EC' = e(elas_c)
forvalues i = 1/4 {
    forvalues j = 1/4 {
        assert abs(_b[ELAS_UNCOMP:e_`i'_`j'] - el(`EU', `i', `j')) < 1e-12
        assert abs(_b[ELAS_COMP:e_`i'_`j'] - el(`EC', `i', `j')) < 1e-12
    }
}
display as result "PASS: every named price elasticity matches its [good, price] matrix cell."

if `test_bootstrap' {
    display as text _newline "=== OPTIONAL BOOTSTRAP SMOKE TEST: reps(4) ==="
    capture noisily pyquaidsce sw1 sw2 sw4 sw9, ///
        prices(p1 p2 p4 p9) expenditure(total) ///
        demographics(x1 x2) anot(10) method(ifgnls) ///
        reps(4) seed(20261011) n_jobs(1) nolog
    local rc_boot = _rc
    if `rc_boot' != 0 {
        display as error "FAIL: bootstrap worker returned Stata code `rc_boot'."
        log close pub_smoke
        exit `rc_boot'
    }
    assert e(N) == 2000
    assert e(converged) == 1
    assert e(boot_reps) >= 2
    display as result "PASS: bootstrap worker returned " e(boot_reps) " successful replications."
}

display as result _newline "SUCCESS: installed Stata smoke test passed."
log close pub_smoke
