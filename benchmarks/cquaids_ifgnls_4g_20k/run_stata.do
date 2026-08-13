* Controlled CQUAIDS/IFGNLS benchmark for Stata quaidsce v2.0
* 20,000 observations, 4 goods, 3 demographics, no initial(), no bootstrap.

clear all
set more off
capture log close _all
set linesize 220

* Run this do-file from benchmarks/cquaids_ifgnls_4g_20k/.
log using "results/stata_benchmark.log", replace text
use "data/benchmark_cquaids_4g_20k.dta", clear

display "STATA CONTROLLED CQUAIDS/IFGNLS BENCHMARK"
display "Stata version: " c(stata_version)
display "OS: " c(os)
display "Machine type: " c(machine_type)
display "Processors available to Stata: " c(processors)
display "N: " _N
foreach v in w1 w2 w3 w4 {
    quietly count if `v' == 0
    display "ZEROS `v' " r(N)
}
display "Specification: censored QUAIDS, IFGNLS, 4 goods, 3 demographics"
display "Starting values: default; no initial() option supplied"
display "Bootstrap: disabled; quaidsce_c is the point-estimation core"

* Data loading is outside the timer. The timer covers the complete point-estimation call.
timer clear 1
timer on 1
capture noisily quaidsce_c w1 w2 w3 w4, anot(10) reps(1)              ///
    prices(p1 p2 p3 p4) expenditure(total) demographics(z1 z2 z3)      ///
    method(ifgnls) nolog
local rc = _rc
timer off 1
quietly timer list 1
scalar bench_runtime = r(t1)
display "BENCHRUNTIME " %21.17g bench_runtime

if `rc' != 0 {
    display "BENCHFAIL return_code=`rc'"
    log close
    exit `rc'
}

* Dump the returned vector and standard errors at full precision.
matrix BB = e(b)
matrix VV = e(V)
local NAMES : colfullnames BB
mata: ///
    b = st_matrix("BB"); v = st_matrix("VV"); nm = tokens(st_local("NAMES")); ///
    for (i=1; i<=cols(b); i++) { ///
        se = (v[i,i] >= 0 ? sqrt(v[i,i]) : .); ///
        printf("BENCHPARM %s %21.17g %21.17g\n", nm[i], b[i], se); ///
    }

display "BENCHLL " %21.17g e(ll)
display "BENCHN " e(N)
log close
