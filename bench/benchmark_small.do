*==============================================================================
* benchmark_small.do  --  a small, fast quaidsce run to validate the Python port
*
* 4 goods, 2000 observations, 2 demographics, censoring in all four shares.
* Should take well under a minute in Stata (no bootstrap: we call quaidsce_c
* directly, so we get the delta-method standard errors, which is exactly what
* the Python port computes).
*
* HOW TO RUN
*   1. put this file next to  DS_STATA_3_2_0_pci2sls_.dta   (it is in the
*      repo's data/ folder, already inside the zip you sent me)
*   2. cd to that folder in Stata and:   do benchmark_small.do
*   3. send me back:   small4.dta  and  small4.log
*
* The mata blocks print every coefficient with 17 significant digits, so I can
* compare against Python at full double precision instead of Stata's 7-digit
* display format.
*==============================================================================

clear all
set more off
capture log close _all
set linesize 200

*------------------------------------------------------------------ build data
use "DS_STATA_3_2_0_pci2sls_.dta", clear

* goods 1, 2, 4 and 9 -- all four have zero shares, so censoring is identified
local gs 1 2 4 9

gen double grp = 0
foreach i of local gs {
    replace grp = grp + w`i'
}
keep if grp > 0                       // conditional (within-group) system
gen double total = total_exp*grp       // group expenditure

foreach i of local gs {
    gen double sw`i' = w`i'/grp        // shares now sum to one
}

keep sw1 sw2 sw4 sw9 p1 p2 p4 p9 total x1 x2
keep if _n <= 2000

summarize
save small4.dta, replace

*---------------------------------------------------------------------- estimate
log using small4.log, replace text

display "==================== RUN 1: default (fgnls) ===================="
quaidsce_c sw1 sw2 sw4 sw9, anot(10) reps(1) prices(p1 p2 p4 p9)      ///
    expenditure(total) demographics(x1 x2) nolog

mat bb = e(b)
mat VV = e(V)
local nm : colfullnames bb
mata: ///
    b = st_matrix("bb"); v = st_matrix("VV"); nmv = tokens(st_local("nm")); ///
    for (i=1; i<=cols(b); i++) ///
        printf("PARM1 %s %21.17g %21.17g\n", nmv[i], b[i], sqrt(v[i,i]))

display "LL1 " %21.17g e(ll)
display "N1 "  e(N)

display "==================== RUN 2: method(nls) ===================="
quaidsce_c sw1 sw2 sw4 sw9, anot(10) reps(1) prices(p1 p2 p4 p9)      ///
    expenditure(total) demographics(x1 x2) nolog method(nls)

mat bb = e(b)
mat VV = e(V)
local nm : colfullnames bb
mata: ///
    b = st_matrix("bb"); v = st_matrix("VV"); nmv = tokens(st_local("nm")); ///
    for (i=1; i<=cols(b); i++) ///
        printf("PARM2 %s %21.17g %21.17g\n", nmv[i], b[i], sqrt(v[i,i]))

display "LL2 " %21.17g e(ll)

display "==================== RUN 3: method(ifgnls) ===================="
quaidsce_c sw1 sw2 sw4 sw9, anot(10) reps(1) prices(p1 p2 p4 p9)      ///
    expenditure(total) demographics(x1 x2) nolog method(ifgnls)

mat bb = e(b)
mat VV = e(V)
local nm : colfullnames bb
mata: ///
    b = st_matrix("bb"); v = st_matrix("VV"); nmv = tokens(st_local("nm")); ///
    for (i=1; i<=cols(b); i++) ///
        printf("PARM3 %s %21.17g %21.17g\n", nmv[i], b[i], sqrt(v[i,i]))

display "LL3 " %21.17g e(ll)

display "==================== RUN 4: noquadratic (AIDS) ===================="
* NOTE: if this run errors out, that is expected and useful information --
* please leave the error in the log and carry on.
capture noisily quaidsce_c sw1 sw2 sw4 sw9, anot(10) reps(1)           ///
    prices(p1 p2 p4 p9) expenditure(total) demographics(x1 x2) nolog noquadratic
if _rc == 0 {
    mat bb = e(b)
    local nm : colfullnames bb
    mata: ///
        b = st_matrix("bb"); nmv = tokens(st_local("nm")); ///
        for (i=1; i<=cols(b); i++) printf("PARM4 %s %21.17g\n", nmv[i], b[i])
    display "LL4 " %21.17g e(ll)
}

display "==================== RUN 5: nocensor (= Poi's quaids) ============"
quaidsce_c sw1 sw2 sw4 sw9, anot(10) reps(1) prices(p1 p2 p4 p9)      ///
    expenditure(total) demographics(x1 x2) nolog nocensor

mat bb = e(b)
mat VV = e(V)
local nm : colfullnames bb
mata: ///
    b = st_matrix("bb"); v = st_matrix("VV"); nmv = tokens(st_local("nm")); ///
    for (i=1; i<=cols(b); i++) ///
        printf("PARM5 %s %21.17g %21.17g\n", nmv[i], b[i], sqrt(v[i,i]))

display "LL5 " %21.17g e(ll)

*--------------------------------------------- first-stage probits, for the record
display "==================== first-stage probits ===================="
foreach i of local gs {
    gen byte z`i' = sw`i' > 0
}
foreach v of varlist p1 p2 p4 p9 total {
    gen double ln`v' = ln(`v')
}
foreach i of local gs {
    display "PROBIT sw`i'"
    probit z`i' lnp1 lnp2 lnp4 lnp9 lntotal x1 x2, nolog
    mat bb = e(b)
    local nm : colfullnames bb
    mata: ///
        b = st_matrix("bb"); nmv = tokens(st_local("nm")); ///
        for (i=1; i<=cols(b); i++) printf("PROBITB %s %21.17g\n", nmv[i], b[i])
}

log close
display "DONE -- please send back small4.dta and small4.log"
