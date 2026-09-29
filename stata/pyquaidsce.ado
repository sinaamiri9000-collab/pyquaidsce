*! version 1.6.0  29sep2026
*! pyquaidsce: backward-compatible alias for the quaidsce Stata command

program define pyquaidsce, eclass
    version 16.0
    quaidsce `0'
end
