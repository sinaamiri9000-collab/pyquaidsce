* pyquaidsce vs quaidsce: 6-food-group comparison
* Run from the root directory of the pyquaidsce repository.
* Install pyquaidsce and quaidsce_c before running.

version 16.0
clear all
set more off
set linesize 120

local DATA "benchmark/pyquaidsce vs quaidsce/six_commodity_input.dta"
local LOG "benchmark/pyquaidsce vs quaidsce/benchmark_new.log"
capture log close comparison
log using "`LOG'", text replace name(comparison)

capture confirm file "`DATA'"
if _rc {
    display as error "Data file not found: `DATA'"
    log close comparison
    exit 601
}
capture noisily which pyquaidsce
if _rc {
    display as error "pyquaidsce is not installed."
    log close comparison
    exit 111
}
capture noisily which quaidsce_c
if _rc {
    display as error "quaidsce_c is not installed."
    log close comparison
    exit 111
}

python query
use "`DATA'", clear
assert _N == 37656
local shares "share1 share2 share3 share4 share5 share6"
local prices "price1 price2 price3 price4 price5 price6"
local demos "age scale children region_1 region_2 region_3 region_4 region_5"
foreach var of varlist `shares' `prices' total_expenditure `demos' {
    assert !missing(`var')
}
foreach var of varlist `prices' total_expenditure {
    assert `var' > 0
}

display as text "=== pyquaidsce ==="
timer clear 1
timer on 1
capture noisily pyquaidsce `shares', prices(`prices') ///
    expenditure(total_expenditure) demographics(`demos') ///
    anot(11) method(ifgnls) algorithm(gn) start(zero) reps(0) nolog ///
    param_tol(1e-5) objective_tol(1e-7) gn_tol(1e-5) ///
    outer_param_tol(1e-5) max_outer(200) max_iter(300) ///
    chunk(2000) blas_threads(1)
local py_rc = _rc
timer off 1
quietly timer list 1
local py_seconds = r(t1)
local py_N = .
local py_ll = .
local py_converged = .
if `py_rc' == 0 {
    local py_N = e(N)
    local py_ll = e(ll)
    local py_converged = e(converged)
}

use "`DATA'", clear
display as text "=== quaidsce_c ==="
timer clear 2
timer on 2
capture noisily quaidsce_c `shares', prices(`prices') ///
    expenditure(total_expenditure) demographics(`demos') ///
    anot(11) method(ifgnls) reps(0) nolog
local st_rc = _rc
timer off 2
quietly timer list 2
local st_seconds = r(t2)
local st_N = .
local st_ll = .
local st_converged = .
if `st_rc' == 0 {
    local st_N = e(N)
    local st_ll = e(ll)
    capture local st_converged = e(converged)
}

display as text "=== Results ==="
display as result "pyquaidsce: rc=`py_rc' N=`py_N' converged=`py_converged' LL=" %17.8f `py_ll' " seconds=" %12.3f `py_seconds'
display as result "quaidsce_c: rc=`st_rc' N=`st_N' converged=`st_converged' LL=" %17.8f `st_ll' " seconds=" %12.3f `st_seconds'
if `py_rc' == 0 & `st_rc' == 0 & `py_N' == `st_N' & `py_seconds' > 0 {
    display as result "Stata / Python time ratio: " %10.3f (`st_seconds' / `py_seconds')
    display as result "Log-likelihood difference: " %12.9g (`py_ll' - `st_ll')
}
else display as error "At least one command failed. Check the errors above."
if `st_converged' >= . & `st_rc' == 0 {
    display as text "Note: quaidsce_c did not report a convergence flag."
}
log close comparison
