* Controlled pyquaidsce / quaidsce runtime benchmark
clear all
set more off
capture log close _all
set linesize 200

* Run this .do file from bench/runtime_comparison
use "benchmark_11goods_37645.dta", clear

assert _N == 37645

log using "stata_runtime.log", replace text

display "STATA_VERSION " c(stata_version)
display "OS " c(os)
display "MACHINE_TYPE " c(machine_type)
display "N " _N

* Time only the estimator. Data loading is deliberately outside the timer.
timer clear 1
timer on 1
quietly quaidsce_c w1 w2 w3 w4 w5 w6 w7 w8 w9 w10 w11, ///
    anot(1.6) reps(1) ///
    prices(tornqvistssb2 tornqvistsweetsnack tornqvistsweetmeal tornqvisttea ///
           tornqvistsoursnack tornqvistfruitveg tornqvistcereals ///
           tornqvistprotein2 tornqvistdairy tornqvistoils tornqvistspices) ///
    expenditure(tfexp) demographics(scale age cfunc) method(ifgnls) nolog
timer off 1

display "FINAL_N " e(N)
display "FINAL_LLF " %21.17g e(ll)
timer list 1

log close
