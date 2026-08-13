*==============================================================================
* validate_my_model.do
*
* Reproduce YOUR paper's model in Python, exactly.
*
* This runs quaidsce_c ONCE on your own data with method(ifgnls) -- no
* bootstrap, so it costs one estimation, not reps() of them -- and writes out
* everything I need to check the Python port against the numbers you will
* actually publish.
*
* EDIT THE FOUR LOCALS BELOW, then:  do validate_my_model.do
* Send me back:  mymodel.dta  and  mymodel.log
*
* mymodel.dta is the *exact estimation sample* after your own cleaning, so
* there is no ambiguity at all about which observations went in.
*==============================================================================

clear all
set more off
capture log close _all
set linesize 250

*------------------------------------------------------------------ EDIT THESE
* your data file (already cleaned, exactly as you feed it to quaidsce)
local DATA      "mydata.dta"

* the variables, in the same order you use them
local SHARES    "w1 w2 w3 w4 w5 w6 w7 w8 w9 w10 w11"
local PRICES    "p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11"
local EXPVAR    "total"
local DEMOS     "x1 x2 x3"
local ANOT      10
*------------------------------------------------------------------------------

use "`DATA'", clear

* keep only what the model needs, and only complete cases -- this is what
* -marksample-/-markout- do inside quaidsce_c, done explicitly so that the
* saved file IS the estimation sample
keep `SHARES' `PRICES' `EXPVAR' `DEMOS'
foreach v of varlist `SHARES' `PRICES' `EXPVAR' `DEMOS' {
    drop if missing(`v')
}
foreach v of varlist `PRICES' `EXPVAR' {
    drop if `v' <= 0
}

* --- do your shares sum to one? -------------------------------------------- *
tempvar swsum
egen double `swsum' = rowtotal(`SHARES')
summarize `swsum'
display "SHARESUM_MEAN " %21.17g r(mean)
display "SHARESUM_MIN  " %21.17g r(min)
display "SHARESUM_MAX  " %21.17g r(max)
* If the mean is not ~1, your system is conditional and the shares should be
* rescaled (w_i/sum_j w_j, with expenditure rescaled to match).  That change
* dramatically improves the numerical conditioning -- see the note I sent.

count
display "NOBS " r(N)
save mymodel.dta, replace

*------------------------------------------------------------------- estimate
log using mymodel.log, replace text

timer clear 1
timer on 1
quaidsce_c `SHARES', anot(`ANOT') reps(1) prices(`PRICES')                ///
    expenditure(`EXPVAR') demographics(`DEMOS') nolog method(ifgnls)
timer off 1
timer list 1

* full double precision dump: 17 significant digits, coefficient and s.e.
mat bb = e(b)
mat VV = e(V)
local nm : colfullnames bb
mata: ///
    b = st_matrix("bb"); v = st_matrix("VV"); nmv = tokens(st_local("nm")); ///
    for (i=1; i<=cols(b); i++) ///
        printf("PARM %s %21.17g %21.17g\n", nmv[i], b[i], sqrt(v[i,i]))

display "LL "     %21.17g e(ll)
display "NOBS "   e(N)
display "ANOT "   %21.17g e(anot)
display "NGOODS " e(ngoods)
display "NDEMOS " e(ndemos)

log close
display "DONE -- send me mymodel.dta and mymodel.log"
