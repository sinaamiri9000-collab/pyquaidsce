* pyquaidsce publication branch: clean installation and small4 smoke tests
* Tested design: Stata 16+ / Windows or macOS/Linux / selected Python 3.9+
* IMPORTANT: Restart Stata before running, so previously imported Python
* modules cannot mask the newly installed publication branch.

version 16.0
clear all
set more off
capture log close publication_check
log using "pyquaidsce_publication_small4_check.log", replace text name(publication_check)

local root "https://raw.githubusercontent.com/sinaamiri9000-collab/pyquaidsce/publication"
local pyexe = c(python_exec)

display as text _newline "=== STATA PYTHON CONFIGURATION ==="
if trim(`"`pyexe'"') == "" {
    display as error "No Python selected. Use: python set exec ...  then restart Stata."
    log close publication_check
    exit 198
}
capture confirm file "`pyexe'"
if _rc {
    display as error "Stata's configured Python executable does not exist: `pyexe'"
    log close publication_check
    exit 601
}
display as text "Python selected by Stata: `pyexe'"
python query

display as text _newline "=== INSTALL PYTHON BACKEND FROM publication ==="
* No git command is required: pip downloads the GitHub branch ZIP.
* --no-deps preserves the user's existing numeric Python environment.
python:
import os
import sys
import subprocess
from sfi import Macro
if "pyquaidsce" in sys.modules:
    raise RuntimeError("pyquaidsce is already imported. Restart Stata and rerun the do-file.")
exe = Macro.getGlobal("c(python_exec)").strip().strip('"')
if not exe or not os.path.isfile(exe):
    raise RuntimeError("No valid Python executable configured in Stata.")
url = "https://github.com/sinaamiri9000-collab/pyquaidsce/archive/refs/heads/publication.zip"
subprocess.run(
    [exe, "-m", "pip", "install", "--upgrade", "--force-reinstall",
     "--no-deps", "--no-cache-dir", url],
    check=True,
)
end

display as text _newline "=== INSTALL STATA FRONTEND FROM publication ==="
net install pyquaidsce, from("https://raw.githubusercontent.com/sinaamiri9000-collab/pyquaidsce/publication/stata") replace
discard
which pyquaidsce

display as text _newline "=== VERIFY PYTHON EXECUTABLE AND SOURCE IDENTITY ==="
python:
import pyquaidsce
from pyquaidsce.stata_bridge import _get_stata_python_executable, _verify_child_python
exe = _get_stata_python_executable()
_verify_child_python(exe)
print("Stata-selected Python:", exe)
print("pyquaidsce version:", pyquaidsce.__version__)
print("pyquaidsce location:", pyquaidsce.__file__)
print("PASS: child Python loads the same package as Stata")
end

display as text _newline "=== DOWNLOAD small4.dta FROM publication ==="
tempfile small4
copy "`root'/bench/small4.dta" "`small4'", replace
use "`small4'", clear

foreach v in sw1 sw2 sw4 sw9 p1 p2 p4 p9 total x1 x2 {
    confirm variable `v'
}
assert _N == 2000
foreach v in p1 p2 p4 p9 total {
    assert `v' > 0 if !missing(`v')
}

* Numeric 0/1 dummy: proves demographic binary regressors remain supported.
generate byte dummy01 = mod(_n, 2)
assert inlist(dummy01, 0, 1)
quietly count if dummy01 == 0
assert r(N) > 0
quietly count if dummy01 == 1
assert r(N) > 0

display as text _newline "=== POINT ESTIMATION / NO BOOTSTRAP ==="
timer clear 1
timer on 1
capture noisily pyquaidsce sw1 sw2 sw4 sw9, ///
    prices(p1 p2 p4 p9) expenditure(total) ///
    demographics(x1 x2 dummy01) anot(10) method(ifgnls) ///
    reps(0) nolog
local rc_point = _rc
timer off 1
if `rc_point' != 0 {
    display as error "FAIL: point estimation (return code `rc_point')."
    log close publication_check
    exit `rc_point'
}
assert e(N) == 2000
assert e(converged) == 1
scalar point_ll = e(ll)
display as result "PASS: point estimation. LL = " %12.5f point_ll
timer list 1

display as text _newline "=== BOOTSTRAP SMOKE TEST (4 replications, NOT inference) ==="
timer clear 2
timer on 2
capture noisily pyquaidsce sw1 sw2 sw4 sw9, ///
    prices(p1 p2 p4 p9) expenditure(total) ///
    demographics(x1 x2 dummy01) anot(10) method(ifgnls) ///
    reps(4) seed(20261011) n_jobs(1) nolog
local rc_boot = _rc
timer off 2
if `rc_boot' != 0 {
    display as error "FAIL: bootstrap run (return code `rc_boot')."
    log close publication_check
    exit `rc_boot'
}
assert e(N) == 2000
assert e(converged) == 1
assert e(boot_reps) >= 2
assert abs(e(ll) - point_ll) < 1e-4
display as result "PASS: bootstrap completed; successful draws = " e(boot_reps) " / 4"
timer list 2

display as result _newline "SUCCESS: publication Python/Stata installation and both smoke tests passed."
display as text "The four bootstrap draws are only a functional test, not valid inference."
log close publication_check
