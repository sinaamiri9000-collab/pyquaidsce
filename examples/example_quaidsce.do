* Make sure the Stata interface and the Python package are installed.
* See stata/README.md for installation instructions.
* Run from the repository root with: do examples/example_quaidsce.do

clear all
set more off

* Read the dataset from your checkout, or download it from the repository.
local data_file "examples/data/household_demand.csv"
capture confirm file "`data_file'"
if _rc local data_file "https://raw.githubusercontent.com/sinaamiri9000-collab/pyquaidsce/main/examples/data/household_demand.csv"
import delimited using "`data_file'", clear

* 1. Basic run: no bootstrap and no expenditure endogeneity correction.
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(expenditure) ///
    demographics(hh_size urban) anot(8) reps(0)
matrix list e(elas_i)
matrix list e(elas_u)
matrix list e(elas_c)

* 2. Use income as an instrument to account for endogenous expenditure.
* The package estimates the expenditure equation and uses its residual.
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(expenditure) ///
    demographics(hh_size urban) anot(8) reps(0) ivexp(income)
matrix list e(reduced_form_b)

* Save starting values from this model for example 5.
matrix b_init = e(b_est)
matrix sigma_init = e(Sigma)

* 3. Run bootstrap with the expenditure correction.
* reps(0) skips bootstrap. reps(10) attempts 10 replications.
* reps(1) attempts one, but at least two must succeed for standard errors.
* We use 10 here for practice; increase reps for your research results.
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(expenditure) ///
    demographics(hh_size urban) anot(8) ivexp(income) reps(10) seed(12345)

* 4. Run bootstrap on two workers at the same time.
* warm starts each replication from the full-sample estimates.
* Each replication may run for at most 900 seconds. A longer replication
* counts as unsuccessful, and the program moves on to the next one.
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(expenditure) ///
    demographics(hh_size urban) anot(8) ivexp(income) reps(10) n_jobs(2) ///
    seed(12345) bootstrap_start(warm) rep_timeout(900)

* 5. Reuse starting values from the same model fitted in example 2.
* Use e(b_est), not e(b): it contains only the free demand parameters.
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(expenditure) ///
    demographics(hh_size urban) anot(8) ivexp(income) reps(0) ///
    initial(b_init) sigma_initial(sigma_init)

* 6. Let the package calculate starting values from a simpler AIDS fit.
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(expenditure) ///
    demographics(hh_size urban) anot(8) start(linear) reps(0)

* 7. Fit AIDS by removing the quadratic expenditure term.
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(expenditure) ///
    demographics(hh_size urban) anot(8) noquadratic reps(0)
