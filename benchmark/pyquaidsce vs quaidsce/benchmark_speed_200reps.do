* ==============================================================================
* BENCHMARK: pyquaidsce vs. quaidsce (C-engine)
* 4-Commodity Censored Quadratic Almost Ideal Demand System (QUAIDS)
* 200 Nonparametric Bootstrap Replications on Parallel CPU Cores
* ==============================================================================
* Description:
*   This do-file executes an end-to-end performance and consistency benchmark
*   comparing the Python-accelerated pyquaidsce package against Stata's native
*   quaidsce command (using the parallel package for multiprocessing).
*   Both estimators are evaluated on the exact same dataset within the same
*   Stata environment to ensure fair, unbiased execution timing.
* ==============================================================================

version 16.0
clear all
set more off

* ------------------------------------------------------------------------------
* 0. Configuration Setup
* ------------------------------------------------------------------------------
local data_file  "food_demand_benchmark_4goods.dta"
local stata_reps "stata_bootstrap_reps_200.dta"

local n_reps     = 200
local n_cores    = 8
local anot_val   = 11.0
local seed_val   = 12345

display _n
display as text "=============================================================================="
display as text "                QUAIDS ESTIMATION BENCHMARK: 200 REPLICATIONS                 "
display as text "=============================================================================="
display as text "Date and Time         : " c(current_date) " " c(current_time)
display as text "Stata Version         : " c(stata_version) " " c(flavor)
display as text "Dataset               : `data_file'"
display as text "Bootstrap Replications: `n_reps'"
display as text "Parallel CPU Cores    : `n_cores'"
display as text "Alpha_0 (anot)        : `anot_val'"
display as text "Random Seed           : `seed_val'"
display as text "==============================================================================" _n

* ------------------------------------------------------------------------------
* 1. Load and Verify Estimation Data
* ------------------------------------------------------------------------------
use "`data_file'", clear

quietly count
local n_obs = r(N)
display as text "Observations in estimation sample: " as result `n_obs'
describe w1 w2 w3 w4 p1 p2 p3 p4 expenditure age scale

* ------------------------------------------------------------------------------
* 2. Benchmark 1: pyquaidsce (Inside Stata)
* ------------------------------------------------------------------------------
display _n
display as text "------------------------------------------------------------------------------"
display as text "BENCHMARK 1: pyquaidsce (`n_reps' bootstrap replications, `n_cores' jobs)"
display as text "------------------------------------------------------------------------------"

timer clear 1
timer on 1

pyquaidsce w1 w2 w3 w4, ///
    prices(p1 p2 p3 p4) expenditure(expenditure) demographics(age scale) ///
    anot(`anot_val') reps(`n_reps') n_jobs(`n_cores') seed(`seed_val')

timer off 1
timer list 1
local t_pyquaidsce = r(t1)
local ll_pyquaidsce = e(ll)
local reps_pyquaidsce = e(boot_reps)

matrix b_py = e(b)
matrix V_py = e(V)

display _n as text ">>> pyquaidsce execution finished in " as result %8.2f `t_pyquaidsce' as text " seconds (" as result %6.2f `t_pyquaidsce'/60 as text " minutes)."

* ------------------------------------------------------------------------------
* 3. Benchmark 2: quaidsce with parallel bs (Stata Native Engine)
* ------------------------------------------------------------------------------
display _n
display as text "------------------------------------------------------------------------------"
display as text "BENCHMARK 2: Stata quaidsce (`n_reps' bootstrap replications, `n_cores' workers)"
display as text "------------------------------------------------------------------------------"

parallel initialize `n_cores'

timer clear 2
timer on 2

parallel bs, reps(`n_reps') saving("`stata_reps'", replace): ///
    quaidsce_c w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(expenditure) ///
    demographics(age scale) anot(`anot_val') reps(0) method(ifgnls) nolog

timer off 2
timer list 2
local t_stata = r(t2)
local ll_stata = e(ll)
local reps_stata = e(N_reps)

matrix b_stata = e(b)
matrix V_stata = e(V)

display _n as text ">>> Stata quaidsce execution finished in " as result %8.2f `t_stata' as text " seconds (" as result %6.2f `t_stata'/60 as text " minutes)."

* ------------------------------------------------------------------------------
* 4. Summary and Performance Comparison
* ------------------------------------------------------------------------------
local speedup = `t_stata' / `t_pyquaidsce'

display _n
display as text "=============================================================================="
display as text "                          BENCHMARK SUMMARY RESULTS                           "
display as text "=============================================================================="
display as text "Metric                            |  pyquaidsce         |  Stata quaidsce     "
display as text "----------------------------------+---------------------+---------------------"
display as text "Software Engine                   |  Python (NumPy/BLAS)|  Stata (C/Mata)     "
display as text "Sample Size (Households)          | " as result %19.0f `n_obs' as text " | " as result %19.0f `n_obs'
display as text "Number of Commodity Groups        | " as result %19.0f 4       as text " | " as result %19.0f 4
display as text "Number of Parameters              | " as result %19.0f 104     as text " | " as result %19.0f 104
display as text "Bootstrap Replications Completed  | " as result %19.0f `reps_pyquaidsce' as text " | " as result %19.0f `reps_stata'
display as text "Parallel CPU Cores Utilized       | " as result %19.0f `n_cores' as text " | " as result %19.0f `n_cores'
display as text "Log-Likelihood                    | " as result %19.4f `ll_pyquaidsce' as text " | " as result %19.4f `ll_stata'
display as text "Execution Time (Seconds)          | " as result %19.2f `t_pyquaidsce' as text " | " as result %19.2f `t_stata'
display as text "Execution Time (Minutes)          | " as result %19.2f `t_pyquaidsce'/60 as text " | " as result %19.2f `t_stata'/60
display as text "----------------------------------+---------------------+---------------------"
display as text "SPEEDUP FACTOR (Stata / Python)   | " as result %17.1f `speedup' as text "x faster" as text " |  1.0x (Baseline)    "
display as text "==============================================================================" _n

* Verify saved Stata replication dataset
display as text "Verifying saved bootstrap replications dataset (`stata_reps')..."
preserve
quietly use "`stata_reps'", clear
quietly describe
display as text "Replication file verified: " as result r(N) as text " replications, " as result r(k) as text " variables recorded."
restore

display as text "Benchmark run completed successfully."

exit, clear
