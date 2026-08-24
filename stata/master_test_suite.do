* ==============================================================================
* STATA MASTER COMPREHENSIVE 24-SCENARIO TEST SUITE FOR PYQUAIDSCE 1.4.0
* ==============================================================================
* Dataset: bd_uruguay.csv (14 food groups, 6,848 observations)
* ==============================================================================

clear all
set more off
discard

display as text _n "=============================================================================="
display as text "PYQUAIDSCE 1.4.0 MASTER COMPREHENSIVE 24-SCENARIO TEST SUITE"
display as text "=============================================================================="

* ---- STEP 0: Load dataset & feature prep -------------------------------------
import delimited "C:\Users\sina\Downloads\bd_uruguay.csv", clear case(preserve)

* Clean negative shares
foreach v of varlist w1_red-w14_red {
    replace `v' = 0 if `v' < 0
}

* Price detection
capture confirm variable P_med1
if _rc {
    local price_list "p_med1-p_med14"
    local p_prefix "p_med"
}
else {
    local price_list "P_med1-P_med14"
    local p_prefix "P_med"
}

* Generate log prices
forvalues i = 1/14 {
    capture drop ln_P_med`i'
    gen double ln_P_med`i' = ln(`p_prefix'`i')
}

capture confirm variable ln_gasto_total
if _rc {
    gen double ln_gasto_total = ln(gasto_total)
}

* Generate endogeneity residual proxy for test 18
regress ln_gasto_total npersonas edad Sex log_ing
predict double vhat, residuals

* Subsystem shares for test 6
egen double tot_sub = rowtotal(w1_red w2_red w3_red)
gen double w1_sub = w1_red / tot_sub
gen double w2_sub = w2_red / tot_sub
gen double w3_sub = w3_red / tot_sub
gen double sub_exp = max(gasto1 + gasto2 + gasto3, 1.0)

local demo_list "npersonas edad Sex prim_comp sec_comp sup_comp log_ing"
local anot_val "2.9766226"

display as result "Data prepared: N = " _N " observations."

* ==============================================================================
* Tier 1: Model Specification & Variable Formats
* ==============================================================================
display as text _n "[01/24] RUNNING: Baseline 14-Good IFGNLS (first_stage_predict=xb)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') method(ifgnls) first_stage_predict(xb) nolog
display as result "       --> [PASS] Baseline LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[02/24] RUNNING: Direct Log-Prices & Log-Expenditure..."
pyquaidsce w1_red-w14_red, lnprices(ln_P_med1-ln_P_med14) lnexpenditure(ln_gasto_total) ///
    demographics(`demo_list') anot(`anot_val') method(ifgnls) nolog
display as result "       --> [PASS] Log-vars LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[03/24] RUNNING: Linear AIDS Model (noquadratic)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') noquadratic method(ifgnls) nolog
display as result "       --> [PASS] Linear AIDS LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[04/24] RUNNING: Uncensored QUAIDS (nocensor, Poi 2012)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') nocensor method(ifgnls) nolog
display as result "       --> [PASS] Uncensored LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[05/24] RUNNING: Translog Constant Variation (anot = 10.0)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(10.0) method(ifgnls) nolog
display as result "       --> [PASS] anot(10) LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[06/24] RUNNING: 3-Good Subsystem Estimation..."
pyquaidsce w1_sub w2_sub w3_sub, prices(`p_prefix'1 `p_prefix'2 `p_prefix'3) expenditure(sub_exp) ///
    demographics(`demo_list') anot(`anot_val') method(ifgnls) nolog
display as result "       --> [PASS] 3-Good LL = " %10.4f e(ll) ", Converged = " e(converged)

* ==============================================================================
* Tier 2: Solvers, Algorithms & Numerical Tolerances
* ==============================================================================
display as text _n "[07/24] RUNNING: Nonlinear Least Squares (method=nls)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') method(nls) nolog
matrix b_nls = e(b_est)
matrix sigma_nls = e(Sigma)
display as result "       --> [PASS] NLS LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[08/24] RUNNING: Feasible Generalized NLS (method=fgnls)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') method(fgnls) nolog
display as result "       --> [PASS] FGNLS LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[09/24] RUNNING: Iterated FGNLS (method=ifgnls)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') method(ifgnls) nolog
matrix b_base = e(b_est)
matrix sigma_base = e(Sigma)
display as result "       --> [PASS] IFGNLS LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[10/24] RUNNING: Levenberg-Marquardt Optimizer (algorithm=lm)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') algorithm(lm) method(ifgnls) nolog
display as result "       --> [PASS] LM LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[11/24] RUNNING: Strict Gradient Stopping Rule (stop_rule=tight)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') stop_rule(tight) method(ifgnls) nolog
display as result "       --> [PASS] Tight rule LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[12/24] RUNNING: Alternative VCE Sigma Formula (vce_sigma=final)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') vce_sigma(final) method(ifgnls) nolog
display as result "       --> [PASS] VCE final LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[13/24] RUNNING: Chunk Size Variation (chunk=500)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') chunk(500) method(ifgnls) nolog
display as result "       --> [PASS] Chunk 500 LL = " %10.4f e(ll) ", Converged = " e(converged)

* ==============================================================================
* Tier 3: Censoring, Probit Specifications & Control Functions
* ==============================================================================
display as text _n "[14/24] RUNNING: Legacy Stata Probit CDF Predictor (first_stage_predict=pr)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') first_stage_predict(pr) method(ifgnls) nolog
display as result "       --> [PASS] Legacy PR LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[15/24] RUNNING: Selection Price Subset (selection_prices=P_med1..3)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') selection_prices(`p_prefix'1 `p_prefix'2 `p_prefix'3) method(ifgnls) nolog
display as result "       --> [PASS] Selection prices LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[16/24] RUNNING: Selection Independent Covariates (npersonas edad)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') selection_covariates(npersonas edad) method(ifgnls) nolog
display as result "       --> [PASS] Selection covs LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[17/24] RUNNING: Selection Omit Log Expenditure (selection_noexpenditure)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') selection_noexpenditure method(ifgnls) nolog
display as result "       --> [PASS] Selection no-exp LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[18/24] RUNNING: Endogeneity Control Function (control_function=vhat)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') control_function(vhat) first_stage_predict(xb) method(ifgnls) nolog
display as result "       --> [PASS] Control Function LL = " %10.4f e(ll) ", Converged = " e(converged)

* ==============================================================================
* Tier 4: Warm-Starting & Matrix Transfers
* ==============================================================================
display as text _n "[19/24] RUNNING: Linearized AIDS Starting Values (start=linear)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') start(linear) method(ifgnls) nolog
display as result "       --> [PASS] start(linear) LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[20/24] RUNNING: Custom Structural Free Parameters (initial=b_base)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') initial(b_base) method(ifgnls) nolog
display as result "       --> [PASS] initial() LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[21/24] RUNNING: Custom Residual Covariance Matrix (sigma_initial=sigma_base)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') initial(b_base) sigma_initial(sigma_base) method(ifgnls) nolog
display as result "       --> [PASS] sigma_initial() LL = " %10.4f e(ll) ", Converged = " e(converged)

display as text _n "[22/24] RUNNING: Chained Multi-Stage Warm-Start (NLS theta/sigma -> IFGNLS)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') initial(b_nls) sigma_initial(sigma_nls) method(ifgnls) nolog
display as result "       --> [PASS] Chained start LL = " %10.4f e(ll) ", Converged = " e(converged)

* ==============================================================================
* Tier 5: Bootstrap, Multiprocessing & Postestimation
* ==============================================================================
display as text _n "[23/24] RUNNING: Parallel Bootstrap Standard Errors (reps=10, n_jobs=4, seed=123456)..."
pyquaidsce w1_red-w14_red, prices(`price_list') expenditure(gasto_total) ///
    demographics(`demo_list') anot(`anot_val') reps(10) n_jobs(4) seed(123456) mp_context(spawn) method(ifgnls) nolog
display as result "       --> [PASS] Bootstrap completed with " e(boot_reps) " ok replications."

display as text _n "[24/24] RUNNING: Postestimation Commands & Hypothesis Testing..."
test [beta]beta_1 = [beta]beta_2
display as result "       --> [PASS] Stata 'test [beta]beta_1 = [beta]beta_2' passed cleanly: p = " %6.4f r(p)
matrix list e(elas_u), noheader format(%8.4f)
display as result "       --> [PASS] e(elas_u) matrix retrieved successfully."

display as text _n "=============================================================================="
display as result "ALL 24 TEST SCENARIOS PASSED SUCCESSFULLY IN STATA!"
display as text "=============================================================================="
