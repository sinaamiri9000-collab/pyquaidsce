*! version 1.7.0  06oct2026
*! quaidsce: recommended common-name alias for the pyquaidsce Stata command

program define quaidsce, eclass
    version 16.0
    pyquaidsce `0'
    ereturn local cmd "quaidsce"
    ereturn local cmdline "quaidsce `0'"
end
