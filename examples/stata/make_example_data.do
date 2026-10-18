*==============================================================================
* make_example_data.do
*
* Generate the reference dataset used by the original Stata example.
*
* The upstream example.do uses Stata's random-number generator. Run this file
* in Stata to obtain the same random-number stream for reference comparisons.
*
* Output: foodcomp_censored_table1.dta
*
* The upstream example.log (27 Nov 2024) differs from example.do in two ways:
*   - the log uses  aux*aux*aux > 20000 , example.do says 10000
*   - the log's -set seed- comes BEFORE -use-
* This file uses the settings recorded in that reference log.
*==============================================================================

clear all
set more off

set seed 123456

use https://www.stata-press.com/data/r18/food_consumption, clear

gen rural  = (runiform() > 0.8)
gen income = exp(rnormal()) + exp(rnormal())*(p_proteins*expfd/10)
gen aux    = (n_adults + n_kids + income - p_fruitveg + p_dairy)*exp(rnormal())

replace w_flours   = 0 if aux*aux*aux > 20000
replace w_proteins = 0 if aux*aux < 30
replace w_fruitveg = 0 if aux < 5
replace w_dairy    = 0 if aux > 30

gen total = w_fruitveg + w_dairy + w_proteins + w_flours
drop aux

* the log reports, in order: 1083 / 1139 / 1035 / 996 changes.
* If your numbers differ, your Stata version's RNG stream differs and the
* comparison against the reference log will not be exact.

* Table 1 of the paper: rescale within the group and drop the outliers, exactly
* as the upstream example.do does, to save the estimation sample.
foreach var in w_proteins w_fruitveg w_dairy w_flours {
    replace `var' = `var'/total
    drop if `var' > .75
    gen d`var' = 0
    replace d`var' = 1 if `var' == 0
}

summarize w* dw*
* example.log reports N = 3,621 with means
*   w_dairy .2183173  w_proteins .3685222  w_fruitveg .2478326  w_flours .1653279

count
save foodcomp_censored_table1.dta, replace
display "Saved foodcomp_censored_table1.dta"
