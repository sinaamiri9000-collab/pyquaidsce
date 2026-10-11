* Simple Stata smoke test for the publication branch.
* First install the Python package OUTSIDE Stata, and install the Stata ado.
* Restart Stata before executing this file.
* No sys.executable overrides, pip calls, or Python path changes are made here.

version 16.0
clear all
set more off

* Set to 1 only AFTER the point-estimate test passes.
local test_bootstrap 0

capture log close pub_smoke
log using "publication_small4_smoke.log", replace text name(pub_smoke)

display as text _newline "=== CHECK PUBLICATION INSTALLATION ==="
capture noisily which pyquaidsce
if _rc {
    display as error "pyquaidsce.ado is not installed. Install it from the publication branch first."
    log close pub_smoke
    exit 111
}

* Version is still 1.7.0, so check for the new publication-only function.
* This import is quick: it does not install anything or estimate a model.
capture noisily python: from pyquaidsce.stata_bridge import _get_stata_python_executable; print("PASS: publication launcher is installed")
if _rc {
    display as error "Publication Python backend not found. Install the publication branch in the Python selected by Stata, then restart Stata."
    log close pub_smoke
    exit 198
}

display as text _newline "=== DOWNLOAD small4.dta ==="
tempfile small4
capture noisily copy "https://raw.githubusercontent.com/sinaamiri9000-collab/pyquaidsce/publication/tests/data/small4.dta" "`small4'", replace
if _rc {
    display as error "Could not download small4.dta; check GitHub access."
    log close pub_smoke
    exit 601
}
use "`small4'", clear
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
display as result "PASS: Stata launched the publication Python worker and received results."
display as result "N = " e(N) "  Log-likelihood = " %16.6f e(ll)
timer list 1

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

display as result _newline "SUCCESS: publication Stata smoke test passed."
log close pub_smoke
