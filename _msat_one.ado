*! version 1.2.0  Araar, Abdelkrim  2026sep08
*! bootstrap helper for msat: re-estimates movestay AND recomputes the effects
*! on each resample, so that every data-dependent step is inside the loop.
*!
*! Usage (mirrors movestay's syntax, plus [xi]):
*!   bootstrap att=r(att) atu=r(atu) ate=r(ate) kappa=r(kappa), reps(200) seed(#): ///
*!       _msat_one depvar indepvars [if] [in], select(treatvar [=] varlist_s) [xi hsize() expand(yes) twostep]

program define _msat_one, rclass
    version 11
    syntax varlist(min=1) [if] [in], SELect(string) [XI HSize(varname) EXPand(string) TWOstep DEPvar(varname)]

    * the treatment variable is the first token of select()
    gettoken treat rest : select, parse(" =")

    if ("`xi'" != "") {
        quietly xi: movestay `varlist' `if' `in', select(`select')
    }
    else {
        quietly movestay `varlist' `if' `in', select(`select')
    }
    local xopts
    if ("`hsize'"  != "") local xopts `xopts' hsize(`hsize')
    if ("`expand'" != "") local xopts `xopts' expand(`expand')
    if ("`twostep'" != "") local xopts `xopts' twostep
    if ("`depvar'" != "") local xopts `xopts' depvar(`depvar')
    quietly msat `treat', `xopts' nose

    return scalar att = r(att)
    return scalar atu = r(atu)
    return scalar ate = r(ate)
    return scalar kappa = r(kappa)
end
