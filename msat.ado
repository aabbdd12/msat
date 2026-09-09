*! version 1.3.2  Araar, Abdelkrim  2026sep09
*! a movestay postcommand: ATT, ATU, ATE, kappa and MTE with delta-method
*! standard errors; expected-effect curves along the quantiles of a ranking
*! variable.
*!
*! Changes from 1.0.0 (2015):
*!   - conditional expectations follow Lokshin & Sajaia (2004) eqs (5)-(8)
*!     [Araar (2015) eqs (10)-(13)]: sign is +lambda1 on the treated
*!     subsample and -lambda0 on the untreated subsample, whichever regime
*!     is predicted. (Fixes the two counterfactual lines of mspredict_ar.)
*!   - standard errors now include parameter-estimation uncertainty via the
*!     delta method on e(V) of movestay.
*!     The 1.0.0 standard errors captured only the sampling variance of the
*!     individual predicted effects and were severely understated.
*!   - r(ate) is now returned (1.0.0 overwrote r(atu) with the ATE).
*!   - no longer depends on mspredict_ar.ado.
*!   - 1.1.1: movestay's e() results are held and restored on exit, so msat
*!     can be called repeatedly (svy: ratio was overwriting them).
*!   - 1.1.2: option nose (point estimates only, for use inside bootstrap).
*!   - 1.2.0: kappa = rho1*sigma1 - rho0*sigma0 with SE (Wald test of
*!            selection on gains); generate() saves the individual expected
*!            effects; rank()/nq() give the expected effect by quantile of a
*!            ranking variable, with SEs and an optional graph; mte() gives
*!            the marginal treatment effect at chosen percentiles of the
*!            participation unobservable, with a common-support flag and an
*!            optional graph.
*!
*!   - 1.3.0: twostep -- re-estimates the model by probit + OLS with the Mills
*!            ratio in each regime (Heckman/Lee two-step) from movestay's
*!            variable lists, and reports the same table; robust to skewed
*!            outcome errors where FIML is not; score, support and MTE from
*!            the probit. semipar -- with mte(), runs a percentile-weights
*!            regression (private engine _msat_pwr, derived from gepwreg,
*!            with the derivative-aware target() bandwidth) ranked on the
*!            probit score and reports the semiparametric MTE next to the
*!            parametric line.
*!   - 1.3.1: the semiparametric results are kept for the graph (return matrix
*!            moved after the graph); engine renamed _msat_pwr to avoid any
*!            collision with the user's gepwreg.
*!   - 1.3.2: supports both movestay naming conventions: the Stata Journal
*!            build (st0071_2, v2.0.0: equations <y>_1 <y>_0, ancillary
*!            lns1 lns2 r1 r2 with index 1 = treated regime, selection
*!            equation named after the treatment variable) and the later
*!            build (<y>0 <y>1, lns0 lns1 r0 r1, equation "select").
*!
*! Syntax:  msat treatvar [if] [in] [, hsize(varname) expand(yes)
*!               nose step(#) level(#) generate(newvar)
*!               rank(varname) nq(#) mte(numlist) pscore(newvar)
*!               twostep semipar depvar(varname)
*!               graph gropts(string)]

program define msat, rclass
    version 11
    syntax varlist(min=1 max=1) [if] [in] [, HSize(varname) EXPand(string) ///
              NOSE STEP(real 0.0001) LEVel(cilevel) GENerate(name)       ///
              RANK(varname) NQ(integer 10) MTE(numlist >0 <1)            ///
              PScore(name) TWOstep SEMIpar DEPvar(varname)               ///
              GRaph GROPTs(string)]
    if ("`e(cmd)'" != "movestay") {
        di as error "msat must follow movestay"
        error 301
    }
    marksample touse
    qui replace `touse' = 0 if !e(sample)
    local treat `varlist'
    if ("`generate'" != "") confirm new variable `generate'
    if ("`pscore'"   != "") confirm new variable `pscore'
    if ("`rank'" == "") local nq = 0
    if (`nq' < 0) {
        di as error "nq() must be positive"
        error 198
    }

    local xopts
    if ("`hsize'"  != "") local xopts `xopts' hsize(`hsize')
    if ("`expand'" != "") local xopts `xopts' expand(`expand')

    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    local p = colsof(`b')
    local depv "`e(depvar)'"
    local wtype "`e(wtype)'"
    local wexp  "`e(wexp)'"
    local mwt
    if ("`wtype'" != "") local mwt "[`wtype'`wexp']"

    * movestay naming conventions (both supported):
    *  (a) Stata Journal build st0071_2 (v2.0.0): e(depvar) = "<treated eq>
    *      <untreated eq> <treatvar>", ancillary lns1 lns2 r1 r2 with index
    *      1 = treated regime, selection equation named after the treatment
    *      variable, regime equations <y>_1 and <y>_0;
    *  (b) later build: e(depvar) = "<untreated eq> <treated eq> select",
    *      ancillary lns0 lns1 r0 r1, regime equations <y>0 and <y>1.
    local cfn : colfullnames e(b)
    local sjnames = (strpos(" `cfn' ", " lns2:_cons ") > 0)
    if (`sjnames') {
        local eq1 : word 1 of `depv'
        local eq0 : word 2 of `depv'
    }
    else {
        local eq0 : word 1 of `depv'
        local eq1 : word 2 of `depv'
    }
    local selname select
    if (strpos(" `cfn' ", " select:") == 0) local selname : word 3 of `depv'
    * e(depvar) re-ordered as "<untreated eq> <treated eq> <selection eq>"
    local depv "`eq0' `eq1' `selname'"

    * variable lists of the three equations, from the column names of e(b)
    local xlist
    local zlist
    foreach nm of local cfn {
        gettoken eqn vn : nm, parse(":")
        local vn = substr("`vn'", 2, .)
        if ("`vn'" == "_cons" | "`vn'" == "") continue
        if ("`eqn'" == "`eq0'")      local xlist `xlist' `vn'
        if ("`eqn'" == "`selname'")  local zlist `zlist' `vn'
    }
    * outcome variable: regime equations are named <depvar>0/<depvar>1 (b)
    * or <depvar>_0/<depvar>_1 (a); a user-specified depvar() overrides
    local yvar "`depvar'"
    if ("`yvar'" == "") {
        if (substr("`eq0'", -2, 2) == "_0") local yvar = substr("`eq0'", 1, length("`eq0'") - 2)
        else                                local yvar = substr("`eq0'", 1, length("`eq0'") - 1)
        capture confirm numeric variable `yvar'
        if (_rc) local yvar
    }
    if ("`twostep'`semipar'" != "" & "`yvar'" == "") {
        di as error "cannot recover the outcome variable from movestay; specify depvar()"
        error 198
    }
    if ("`twostep'" != "") local xopts `xopts' ts

    * svy: ratio below overwrites e(); hold movestay's results and restore
    * them automatically when msat exits (also on error), so that msat can
    * be called repeatedly and other postestimation commands still work.
    tempname eshold
    _estimates hold `eshold', restore copy

    * ---- twostep: probit + OLS with the Mills ratio in each regime --------
    tempvar zg2
    if ("`twostep'" != "") {
        _msat_twostep `treat' if `touse', yvar(`yvar') xlist(`xlist') zlist(`zlist') ///
            eq0(`eq0') eq1(`eq1') wt(`mwt') zg(`zg2')
        matrix `b' = r(b)
        matrix `V' = r(V)
        local p = colsof(`b')
    }

    * ------------------------------------------------------------------
    * 1. individual expected effects, observable gain, propensity score
    * ------------------------------------------------------------------
    tempvar d m ps hs
    qui _msat_eff `treat' if `touse', b(`b') gen(`d') mgen(`m') pgen(`ps') ///
        depvar(`depv') `xopts'
    local kappa = r(kappa)
    qui gen double `hs' = 1
    if ("`hsize'" != "") qui replace `hs' = `hsize'

    * aggregation weight for all weighted sums: survey pweight (from svyset,
    * if any) times the size variable. svy: ratio applies the pweight itself
    * and therefore receives `hs' only.
    qui svyset
    local wvar "`r(wvar)'"
    if ("`r(settings)'" == ", clear") {
        qui svyset _n, vce(linearized)
        local wvar
    }
    tempvar hw
    qui gen double `hw' = `hs'
    if ("`wvar'" != "") qui replace `hw' = `hs' * `wvar'

    if ("`generate'" != "") {
        qui gen double `generate' = `d'
        label variable `generate' "E[Y1-Y0 | D, X, Z] from msat"
    }
    if ("`pscore'" != "") {
        qui gen double `pscore' = `ps'
        if ("`twostep'" != "") label variable `pscore' "P(D=1 | Z) from probit (msat twostep)"
        else                   label variable `pscore' "P(D=1 | Z) from movestay, via msat"
    }

    * common support of the propensity score P(Z) between the two groups
    qui sum `ps' if `touse' & `treat'==1
    local p1lo = r(min)
    local p1hi = r(max)
    qui sum `ps' if `touse' & `treat'==0
    local p0lo = r(min)
    local p0hi = r(max)
    local plo = max(`p1lo', `p0lo')
    local phi = min(`p1hi', `p0hi')

    * quantile bins of the ranking variable (whole estimation sample)
    tempvar q
    if (`nq' > 0) {
        qui xtile `q' = `rank' [aw=`hw'] if `touse', nq(`nq')
        local qopts q(`q') nq(`nq')
    }
    else local qopts

    * all point estimates in one row vector S :
    *   [att atu ate m kappa | att_1 atu_1 all_1 ... att_nq atu_nq all_nq]
    tempname S
    _msat_sums, treat(`treat') touse(`touse') d(`d') hs(`hw') m(`m') ///
                k(`kappa') `qopts'
    matrix `S' = r(S)
    local K = colsof(`S')
    local att = el(`S',1,1)
    local atu = el(`S',1,2)
    local ate = el(`S',1,3)
    local mg  = el(`S',1,4)

    * ---- nose : point estimates only (fast path, e.g. inside bootstrap) ----
    if ("`nose'" != "") {
        if ("`twostep'" != "") di as txt _n "Treatment effects, two-step switching regression (msat 1.3)"
        else                   di as txt _n "Treatment effects after endogenous switching regression (msat 1.3)"
        di as txt "{hline 9}{c TT}{hline 14}"
        di as txt "Index    {c |}   Estimate"
        di as txt "{hline 9}{c +}{hline 14}"
        foreach k in att atu ate {
            di as txt %8s upper("`k'") " {c |}" as res %11.6g ``k''
        }
        di as txt %8s "kappa" " {c |}" as res %11.6g `kappa'
        di as txt "{hline 9}{c BT}{hline 14}"
        return scalar att   = `att'
        return scalar atu   = `atu'
        return scalar ate   = `ate'
        return scalar kappa = `kappa'
        return scalar m     = `mg'
        exit
    }

    * ------------------------------------------------------------------
    * 2. sampling-variance component (design-based)
    * ------------------------------------------------------------------
    tempname VS
    matrix `VS' = J(1, `K', 0)
    qui svy: ratio `d'/`hs' if `touse' & `treat'==1
    matrix `VS'[1,1] = el(e(V),1,1)
    qui svy: ratio `d'/`hs' if `touse' & `treat'==0
    matrix `VS'[1,2] = el(e(V),1,1)
    qui svy: ratio `d'/`hs' if `touse'
    matrix `VS'[1,3] = el(e(V),1,1)
    qui svy: ratio `m'/`hs' if `touse'
    matrix `VS'[1,4] = el(e(V),1,1)
    forvalues j = 1/`nq' {
        local c = 5 + 3*(`j'-1)
        qui count if `touse' & `treat'==1 & `q'==`j'
        if (r(N) > 1) {
            qui svy: ratio `d'/`hs' if `touse' & `treat'==1 & `q'==`j'
            matrix `VS'[1,`c'+1] = el(e(V),1,1)
        }
        qui count if `touse' & `treat'==0 & `q'==`j'
        if (r(N) > 1) {
            qui svy: ratio `d'/`hs' if `touse' & `treat'==0 & `q'==`j'
            matrix `VS'[1,`c'+2] = el(e(V),1,1)
        }
        qui count if `touse' & `q'==`j'
        if (r(N) > 1) {
            qui svy: ratio `d'/`hs' if `touse' & `q'==`j'
            matrix `VS'[1,`c'+3] = el(e(V),1,1)
        }
    }

    * ------------------------------------------------------------------
    * 3. delta-method component: G * V * G'
    * ------------------------------------------------------------------
    tempname G bp bm Sp Sm
    matrix `G' = J(`K', `p', 0)
    forvalues j = 1/`p' {
        local bj = el(`b',1,`j')
        local h  = `step' * max(1, abs(`bj'))
        matrix `bp' = `b'
        matrix `bm' = `b'
        matrix `bp'[1,`j'] = `bj' + `h'
        matrix `bm'[1,`j'] = `bj' - `h'

        tempvar dp dm mp mm
        qui _msat_eff `treat' if `touse', b(`bp') gen(`dp') mgen(`mp') ///
            depvar(`depv') `xopts'
        local kp = r(kappa)
        qui _msat_eff `treat' if `touse', b(`bm') gen(`dm') mgen(`mm') ///
            depvar(`depv') `xopts'
        local km = r(kappa)
        _msat_sums, treat(`treat') touse(`touse') d(`dp') hs(`hw') m(`mp') ///
                    k(`kp') `qopts'
        matrix `Sp' = r(S)
        _msat_sums, treat(`treat') touse(`touse') d(`dm') hs(`hw') m(`mm') ///
                    k(`km') `qopts'
        matrix `Sm' = r(S)
        forvalues i = 1/`K' {
            local gij = (el(`Sp',1,`i') - el(`Sm',1,`i')) / (2*`h')
            if (`gij' < .) matrix `G'[`i',`j'] = `gij'
        }
        drop `dp' `dm' `mp' `mm'
    }
    tempname VD
    matrix `VD' = `G' * `V' * `G''

    * total standard errors
    tempname SE
    matrix `SE' = J(1, `K', .)
    forvalues i = 1/`K' {
        if (el(`S',1,`i') < .) {
            matrix `SE'[1,`i'] = sqrt(el(`VS',1,`i') + el(`VD',`i',`i'))
        }
    }
    local se_att   = el(`SE',1,1)
    local se_atu   = el(`SE',1,2)
    local se_ate   = el(`SE',1,3)
    local se_m     = el(`SE',1,4)
    local se_kappa = el(`SE',1,5)
    local vd_att = el(`VD',1,1)
    local vd_atu = el(`VD',2,2)
    local vd_ate = el(`VD',3,3)
    local vs_att = el(`VS',1,1)
    local vs_atu = el(`VS',1,2)
    local vs_ate = el(`VS',1,3)

    * ------------------------------------------------------------------
    * 4. display: main table
    * ------------------------------------------------------------------
    local z = invnormal(1 - (100-`level')/200)
    if ("`twostep'" != "") {
        di as txt _n "Treatment effects, two-step switching regression (msat 1.3)"
        di as txt "Estimation: probit, then OLS on X and the Mills ratio in each regime"
        di as txt "Std. err.: delta method on the two-step estimates (block-diagonal across"
        di as txt "           stages) + sampling component; bootstrap with _msat_one for exact SEs"
    }
    else {
        di as txt _n "Treatment effects after endogenous switching regression (msat 1.3)"
        di as txt "Std. err.: delta method on e(V) + design-based sampling component"
    }
    di as txt "{hline 9}{c TT}{hline 66}"
    di as txt "Index    {c |}   Estimate   Std. err.        z    P>|z|     [`level'% conf. interval]"
    di as txt "{hline 9}{c +}{hline 66}"
    foreach k in att atu ate kappa {
        local est = ``k''
        local se  = `se_`k''
        local zz  = `est'/`se'
        local pv  = 2*(1-normal(abs(`zz')))
        local lb  = `est' - `z'*`se'
        local ub  = `est' + `z'*`se'
        if ("`k'" == "kappa") local lab "kappa"
        else                  local lab = upper("`k'")
        di as txt %8s "`lab'" " {c |}" as res %11.6g `est' %12.6g `se' %9.2f `zz' ///
           %9.3f `pv' %12.6g `lb' %11.6g `ub'
    }
    di as txt "{hline 9}{c BT}{hline 66}"
    if ("`twostep'" != "") di as txt "kappa = coefficient on the Mills ratio, regime 1 minus regime 0; " ///
              "ATT - ATU = kappa*(mean lambda1 + mean lambda0)"
    else di as txt "kappa = rho1*sigma1 - rho0*sigma0 (selection on gains); " ///
              "ATT - ATU = kappa*(mean lambda1 + mean lambda0)"
    di as txt "Variance components (ATT): parameter " %9.3g `vd_att' ///
              "   sampling " %9.3g `vs_att'
    di as txt "Common support of P(Z): [" %6.4f `plo' ", " %6.4f `phi' "]"

    * ------------------------------------------------------------------
    * 5. expected effect by quantile of the ranking variable
    * ------------------------------------------------------------------
    if (`nq' > 0) {
        tempname CV
        matrix `CV' = J(`nq', 7, .)
        di as txt _n "Expected treatment effect by quantile of `rank' (nq = `nq')"
        di as txt "{hline 5}{c TT}{hline 24}{c TT}{hline 24}{c TT}{hline 24}"
        di as txt "  q  {c |}   Treated    Std. err.{c |}  Untreated   Std. err.{c |}    All       Std. err."
        di as txt "{hline 5}{c +}{hline 24}{c +}{hline 24}{c +}{hline 24}"
        forvalues j = 1/`nq' {
            local c = 5 + 3*(`j'-1)
            matrix `CV'[`j',1] = `j'
            forvalues g = 1/3 {
                local e`g' = el(`S',1,`c'+`g')
                local s`g' = el(`SE',1,`c'+`g')
                matrix `CV'[`j', 2*`g']   = `e`g''
                matrix `CV'[`j', 2*`g'+1] = `s`g''
            }
            di as txt %4.0f `j' " {c |}" as res %11.4f `e1' %12.4f `s1' as txt "{c |}" ///
               as res %11.4f `e2' %12.4f `s2' as txt "{c |}" ///
               as res %11.4f `e3' %12.4f `s3'
        }
        di as txt "{hline 5}{c BT}{hline 24}{c BT}{hline 24}{c BT}{hline 24}"
        matrix colnames `CV' = q att se_att atu se_atu all se_all
        return matrix curve = `CV'
    }

    * ------------------------------------------------------------------
    * 6. marginal treatment effect at chosen percentiles u of the
    *    participation unobservable:  MTE(u) = m + kappa * invnormal(1-u)
    * ------------------------------------------------------------------
    if ("`mte'" != "") {
        local vm  = el(`VS',1,4) + el(`VD',4,4)
        local vk  = el(`VD',5,5)
        local vmk = el(`VD',4,5)
        local nm : word count `mte'
        tempname MT
        matrix `MT' = J(`nm', 4, .)
        di as txt _n "Marginal treatment effect at percentiles of the participation unobservable"
        if ("`expand'" == "yes") di as txt "(reported on the scale of the linear index)"
        di as txt "{hline 8}{c TT}{hline 24}{c TT}{hline 12}"
        di as txt "   u    {c |}   MTE(u)     Std. err.{c |}  support"
        di as txt "{hline 8}{c +}{hline 24}{c +}{hline 12}"
        local i = 0
        foreach u of numlist `mte' {
            local ++i
            local cu  = invnormal(1 - `u')
            local est = `mg' + `kappa'*`cu'
            local se  = sqrt(`vm' + (`cu')^2*`vk' + 2*(`cu')*`vmk')
            local ins = (`u' >= `plo' & `u' <= `phi')
            matrix `MT'[`i',1] = `u'
            matrix `MT'[`i',2] = `est'
            matrix `MT'[`i',3] = `se'
            matrix `MT'[`i',4] = `ins'
            if (`ins') local tag "  observed"
            else       local tag "  extrapol."
            di as txt %7.3f `u' " {c |}" as res %11.4f `est' %12.4f `se' as txt "{c |}`tag'"
        }
        di as txt "{hline 8}{c BT}{hline 24}{c BT}{hline 12}"
        di as txt "support: u inside the common support of P(Z) " ///
                  "[" %6.4f `plo' ", " %6.4f `phi' "]"
        matrix colnames `MT' = u mte se support
        return matrix mte = `MT'
    }

    * ------------------------------------------------------------------
    * 6b. semiparametric MTE by percentile-weights regression on the
    *     probit score (gepwreg 1.4+, target()), at the mte() percentiles
    * ------------------------------------------------------------------
    local nsp = 0
    if ("`semipar'" != "" & "`mte'" != "") {
        capture which _msat_pwr
        if (_rc) {
            di as txt _n "semipar: _msat_pwr.ado not found (it ships with msat) -- semiparametric MTE skipped"
        }
        else {
            tempvar pzs pz2
            if ("`twostep'" != "") qui gen double `pzs' = normal(`zg2') if `touse'
            else {
                qui probit `treat' `zlist' `mwt' if `touse'
                qui predict double `pzs' if `touse', pr
            }
            qui gen double `pz2' = `pzs'^2 if `touse'
            local wsp
            local mwta = subinstr("`mwt'", "pweight", "aweight", 1)
            foreach v of local xlist {
                tempvar c_`v' i_`v'
                qui sum `v' `mwta' if `touse', meanonly
                qui gen double `c_`v'' = `v' - r(mean) if `touse'
                qui gen double `i_`v'' = `c_`v''*`pzs' if `touse'
                local wsp `wsp' `c_`v'' `i_`v''
            }
            local nsp : word count `mte'
            tempname SP
            matrix `SP' = J(`nsp', 8, .)
            local i = 0
            foreach u of numlist `mte' {
                local ++i
                capture qui _msat_pwr `yvar' `wsp' `pzs' `pz2' `mwt' if `touse', ///
                    per(`u') rankvar(`pzs') target(`pzs' `pz2')
                if (_rc) {
                    di as txt "semipar: percentile-weights regression failed at tau = `u' (rc = " _rc ")"
                    continue
                }
                local qt = e(q_tau)
                local hh = e(h)
                local ne = e(N_eff)
                local lo = max(`u' - `hh'*sqrt(2), 0.005)
                local hi = min(`u' + `hh'*sqrt(2), 0.995)
                _pctile `pzs' if `touse', p(`=100*`lo'' `=100*`hi'')
                local wid = r(r2) - r(r1)
                qui lincom _b[`pzs'] + 2*`qt'*_b[`pz2']
                matrix `SP'[`i',1] = `u'
                matrix `SP'[`i',2] = `qt'
                matrix `SP'[`i',3] = `mg' + `kappa'*invnormal(1-`qt')
                matrix `SP'[`i',4] = r(estimate)
                matrix `SP'[`i',5] = r(se)
                matrix `SP'[`i',6] = `hh'
                matrix `SP'[`i',7] = `ne'
                matrix `SP'[`i',8] = `wid'
            }
            di as txt _n "Semiparametric MTE (percentile-weights regression on the probit score, derivative-aware bandwidth)"
            di as txt "{hline 7}{c TT}{hline 8}{c TT}{hline 11}{c TT}{hline 22}{c TT}{hline 8}{c TT}{hline 8}{c TT}{hline 9}"
            di as txt "  tau  {c |} p_tau  {c |} parametric{c |}     PWR     Std. err.{c |}    h   {c |}  N_eff {c |} width p"
            di as txt "{hline 7}{c +}{hline 8}{c +}{hline 11}{c +}{hline 22}{c +}{hline 8}{c +}{hline 8}{c +}{hline 9}"
            forvalues j = 1/`nsp' {
                di as txt %6.2f `SP'[`j',1] " {c |}" as res %7.3f `SP'[`j',2] as txt " {c |}" ///
                   as res %10.4f `SP'[`j',3] as txt " {c |}" ///
                   as res %10.4f `SP'[`j',4] %11.4f `SP'[`j',5] as txt " {c |}" ///
                   as res %7.3f `SP'[`j',6] as txt " {c |}" as res %7.0f `SP'[`j',7] as txt " {c |}" ///
                   as res %8.3f `SP'[`j',8]
            }
            di as txt "{hline 7}{c BT}{hline 8}{c BT}{hline 11}{c BT}{hline 22}{c BT}{hline 8}{c BT}{hline 8}{c BT}{hline 9}"
            di as txt "tau: quantile of the probit score; p_tau: the score at that quantile; parametric: m + kappa*invnormal(1-p_tau)."
            di as txt "PWR std. err. exclude the estimation of the score. A gap between the two columns that grows"
            di as txt "towards the tails is the signature of a departure from joint normality."
            matrix colnames `SP' = tau p_tau parametric pwr se h neff width_p
        }
    }

    * ------------------------------------------------------------------
    * 7. graphs
    * ------------------------------------------------------------------
    if ("`graph'" != "") {
        preserve
        if (`nq' > 0) {
            qui drop _all
            qui set obs `nq'
            qui gen q = _n
            foreach v in att se_att atu se_atu all se_all {
                qui gen double `v' = .
            }
            forvalues j = 1/`nq' {
                local c = 5 + 3*(`j'-1)
                qui replace att    = el(`S',1,`c'+1)  in `j'
                qui replace se_att = el(`SE',1,`c'+1) in `j'
                qui replace atu    = el(`S',1,`c'+2)  in `j'
                qui replace se_atu = el(`SE',1,`c'+2) in `j'
                qui replace all    = el(`S',1,`c'+3)  in `j'
                qui replace se_all = el(`SE',1,`c'+3) in `j'
            }
            foreach v in att atu all {
                qui gen double `v'_lb = `v' - `z'*se_`v'
                qui gen double `v'_ub = `v' + `z'*se_`v'
            }
            twoway (rarea att_lb att_ub q, color(gs13) lwidth(none))          ///
                   (rarea atu_lb atu_ub q, color(gs15) lwidth(none))          ///
                   (line att q, lpattern(solid) lcolor(black))                ///
                   (line atu q, lpattern(dash)  lcolor(black))                ///
                   (line all q, lpattern(solid) lcolor(gs8) lwidth(medthick)) ///
                   , legend(order(3 "Treated (ATT)" 4 "Untreated (ATU)" 5 "All (ATE)") ///
                     rows(1)) ytitle("Expected treatment effect")             ///
                     xtitle("Quantile of `rank'") xlabel(1(1)`nq')            ///
                     title("Expected effect by quantile of `rank'")          ///
                     note("`level'% confidence bands: delta method + sampling") ///
                     name(msat_rank, replace) `gropts'
        }
        if ("`mte'" != "") {
            qui drop _all
            qui set obs 99
            qui gen double u  = _n/100
            qui gen double cu = invnormal(1-u)
            qui gen double mte = `mg' + `kappa'*cu
            qui gen double se  = sqrt(`vm' + cu^2*`vk' + 2*cu*`vmk')
            qui gen double lb  = mte - `z'*se
            qui gen double ub  = mte + `z'*se
            local spplot
            local splegend
            if (`nsp' > 0) {
                qui set obs `=99+`nsp''
                qui gen double qsp  = .
                qui gen double msp  = .
                qui gen double lbsp = .
                qui gen double ubsp = .
                forvalues j = 1/`nsp' {
                    qui replace qsp  = `SP'[`j',2] in `=99+`j''
                    qui replace msp  = `SP'[`j',4] in `=99+`j''
                    qui replace lbsp = `SP'[`j',4] - `z'*`SP'[`j',5] in `=99+`j''
                    qui replace ubsp = `SP'[`j',4] + `z'*`SP'[`j',5] in `=99+`j''
                }
                local spplot (rcap lbsp ubsp qsp, lcolor(gs6)) (scatter msp qsp, mcolor(black) msymbol(O))
                local splegend 4 "PWR on the probit score (semiparametric)"
            }
            twoway (rarea lb ub u, color(gs13) lwidth(none))                  ///
                   (line mte u, lcolor(black))                                ///
                   `spplot'                                                   ///
                   , xline(`plo' `phi', lpattern(shortdash) lcolor(gs8))      ///
                     yline(`ate', lpattern(dot) lcolor(gs8))                  ///
                     legend(order(2 "MTE(u), parametric" `splegend') rows(1)) ///
                     ytitle("Marginal treatment effect")                      ///
                     xtitle("u : percentile of the participation unobservable (0 = most eager)") ///
                     title("Marginal treatment effect")                       ///
                     note("Dashed verticals: common support of P(Z). Dotted: ATE.") ///
                     name(msat_mte, replace) `gropts'
        }
        restore
    }

    * ------------------------------------------------------------------
    * 8. returns
    * ------------------------------------------------------------------
    tempname PS
    matrix `PS' = (`plo', `phi')
    matrix colnames `PS' = lo hi
    return matrix psupport = `PS'
    if (`nsp' > 0) return matrix mte_sp = `SP'
    return matrix G = `G'
    return scalar att      = `att'
    return scalar atu      = `atu'
    return scalar ate      = `ate'
    return scalar kappa    = `kappa'
    return scalar m        = `mg'
    return scalar se_att   = `se_att'
    return scalar se_atu   = `se_atu'
    return scalar se_ate   = `se_ate'
    return scalar se_kappa = `se_kappa'
    return scalar se_m     = `se_m'
    return scalar vd_att = `vd_att'
    return scalar vd_atu = `vd_atu'
    return scalar vd_ate = `vd_ate'
    return scalar vs_att = `vs_att'
    return scalar vs_atu = `vs_atu'
    return scalar vs_ate = `vs_ate'
    return scalar p_lo   = `plo'
    return scalar p_hi   = `phi'
end


* ======================================================================
* _msat_sums : weighted means of the individual effects, overall and by
*              quantile bin, packed in one row vector
*   S = [att atu ate m kappa | att_1 atu_1 all_1 ... att_nq atu_nq all_nq]
* ======================================================================
program define _msat_sums, rclass
    version 11
    syntax , TREAT(varname) TOUSE(varname) D(varname) HS(varname) M(varname) ///
             K(real) [Q(varname) NQ(integer 0)]
    tempname S
    local K = 5 + 3*`nq'
    matrix `S' = J(1, `K', .)
    qui sum `d' [aw=`hs'] if `touse' & `treat'==1
    matrix `S'[1,1] = r(mean)
    qui sum `d' [aw=`hs'] if `touse' & `treat'==0
    matrix `S'[1,2] = r(mean)
    qui sum `d' [aw=`hs'] if `touse'
    matrix `S'[1,3] = r(mean)
    qui sum `m' [aw=`hs'] if `touse'
    matrix `S'[1,4] = r(mean)
    matrix `S'[1,5] = `k'
    forvalues j = 1/`nq' {
        local c = 5 + 3*(`j'-1)
        qui sum `d' [aw=`hs'] if `touse' & `treat'==1 & `q'==`j'
        if (r(N) > 0) matrix `S'[1,`c'+1] = r(mean)
        qui sum `d' [aw=`hs'] if `touse' & `treat'==0 & `q'==`j'
        if (r(N) > 0) matrix `S'[1,`c'+2] = r(mean)
        qui sum `d' [aw=`hs'] if `touse' & `q'==`j'
        if (r(N) > 0) matrix `S'[1,`c'+3] = r(mean)
    }
    return matrix S = `S'
end


* ======================================================================
* _msat_eff : individual treatment effects d_i from a parameter vector b
*             (same layout and names as movestay's e(b))
*   d_i = E[Y1|D_i] - E[Y0|D_i], with the conditioning subsample fixing
*   the sign of the Mills term:  +lambda1 if D=1,  -lambda0 if D=0.
*   Optionally: mgen() = x'(b1-b0) (observable gain), pgen() = P(Z);
*   returns r(kappa) = rho1*sigma1 - rho0*sigma0.
* ======================================================================
program define _msat_eff, rclass
    version 11
    syntax varlist(min=1 max=1) [if] [in], B(name) GEN(name) DEPvar(string) ///
                              [MGen(name) PGen(name) HSize(varname) EXPand(string) TS]
    marksample touse
    local treat `varlist'

    tempname TMP s0 s1 r0 r1 B0 B1 G
    local cn : colfullnames `b'
    if ("`ts'" == "") {
        if (strpos(" `cn' ", " lns2:_cons ") > 0) {
            * Stata Journal build: index 1 = treated, 2 = untreated
            local n_s0 lns2
            local n_s1 lns1
            local n_r0 r2
            local n_r1 r1
        }
        else {
            local n_s0 lns0
            local n_s1 lns1
            local n_r0 r0
            local n_r1 r1
        }
        matrix `TMP' = `b'[1,"`n_s0':_cons"]
        scalar `s0' = exp(`TMP'[1,1])
        matrix `TMP' = `b'[1,"`n_s1':_cons"]
        scalar `s1' = exp(`TMP'[1,1])
        matrix `TMP' = `b'[1,"`n_r0':_cons"]
        scalar `r0' = tanh(`TMP'[1,1])
        matrix `TMP' = `b'[1,"`n_r1':_cons"]
        scalar `r1' = tanh(`TMP'[1,1])
    }
    else {
        * two-step: rs_j = coefficient on the Mills ratio = rho_j*sigma_j
        matrix `TMP' = `b'[1,"rs0:_cons"]
        scalar `r0' = `TMP'[1,1]
        scalar `s0' = 1
        matrix `TMP' = `b'[1,"rs1:_cons"]
        scalar `r1' = `TMP'[1,1]
        scalar `s1' = 1
    }

    * depvar() = "eq0 eq1 seleq" as re-ordered by msat: eq0 = regime 0
    * (untreated), eq1 = regime 1, seleq = name of the selection equation
    * (the two-step vector always names it "select")
    local y0 : word 1 of `depvar'
    local y1 : word 2 of `depvar'
    local ys : word 3 of `depvar'
    if (strpos(" `cn' ", " select:") > 0) local ys select
    matrix `B0' = `b'[1,"`y0':"]
    matrix `B1' = `b'[1,"`y1':"]
    matrix `G'  = `b'[1,"`ys':"]

    tempvar xb0 xb1 zg l1 l0
    qui matrix score double `xb0' = `B0' if `touse'
    qui matrix score double `xb1' = `B1' if `touse'
    qui matrix score double `zg'  = `G'  if `touse'
    qui gen double `l1' = normalden(`zg')/normal(`zg')     if `touse'
    qui gen double `l0' = normalden(`zg')/(1-normal(`zg')) if `touse'

    qui gen double `gen' = .
    if ("`expand'" == "yes") {
        qui replace `gen' = exp(`xb1' + `r1'*`s1'*`l1') - exp(`xb0' + `r0'*`s0'*`l1') ///
            if `touse' & `treat'==1
        qui replace `gen' = exp(`xb1' - `r1'*`s1'*`l0') - exp(`xb0' - `r0'*`s0'*`l0') ///
            if `touse' & `treat'==0
    }
    else {
        qui replace `gen' = (`xb1' + `r1'*`s1'*`l1') - (`xb0' + `r0'*`s0'*`l1') ///
            if `touse' & `treat'==1
        qui replace `gen' = (`xb1' - `r1'*`s1'*`l0') - (`xb0' - `r0'*`s0'*`l0') ///
            if `touse' & `treat'==0
    }
    if ("`mgen'" != "") qui gen double `mgen' = `xb1' - `xb0' if `touse'
    if ("`pgen'" != "") qui gen double `pgen' = normal(`zg')   if `touse'
    return scalar kappa = `r1'*`s1' - `r0'*`s0'
end


* ======================================================================
* _msat_twostep : probit + OLS with the Mills ratio in each regime.
*   Returns r(b), r(V) laid out with movestay's equation names so that
*   _msat_eff (ts mode) and the delta loop can be reused unchanged:
*     [eq0: X _cons | eq1: X _cons | select: Z _cons | rs0:_cons | rs1:_cons]
*   V is block-diagonal across the three estimations (regime 0 with rs0,
*   regime 1 with rs1, probit); it ignores the covariance induced by the
*   generated regressor, hence the bootstrap recommendation for exact SEs.
* ======================================================================
program define _msat_twostep, rclass
    version 11
    syntax varlist(min=1 max=1) [if] [in], YVAR(varname) XLIST(string) ZLIST(string) ///
                              EQ0(string) EQ1(string) ZG(name) [WT(string)]
    marksample touse
    local treat `varlist'
    local kx : word count `xlist'
    local kz : word count `zlist'
    local kx1 = `kx' + 1
    local kz1 = `kz' + 1
    local p = 2*`kx1' + `kz1' + 2

    * probit
    qui probit `treat' `zlist' `wt' if `touse'
    tempname bg Vg
    matrix `bg' = e(b)
    matrix `Vg' = e(V)
    qui predict double `zg' if `touse', xb
    tempvar l1 l0 ml
    qui gen double `l1' = normalden(`zg')/normal(`zg')     if `touse'
    qui gen double `l0' = normalden(`zg')/(1-normal(`zg')) if `touse'
    qui gen double `ml' = cond(`treat'==1, `l1', -`l0')    if `touse'

    * regime regressions (robust): columns X..., ml, _cons
    tempname b1 V1 b0 V0
    qui regress `yvar' `xlist' `ml' `wt' if `touse' & `treat'==1, vce(robust)
    matrix `b1' = e(b)
    matrix `V1' = e(V)
    qui regress `yvar' `xlist' `ml' `wt' if `touse' & `treat'==0, vce(robust)
    matrix `b0' = e(b)
    matrix `V0' = e(V)

    * assemble b and V
    tempname B VV
    matrix `B'  = J(1, `p', 0)
    matrix `VV' = J(`p', `p', 0)
    * index maps: regress order (1..kx = X, kx+1 = ml, kx+2 = _cons) -> joint order
    *   regime 0: X_j -> j ; _cons -> kx1 ; ml -> 2*kx1 + kz1 + 1
    *   regime 1: X_j -> kx1 + j ; _cons -> 2*kx1 ; ml -> 2*kx1 + kz1 + 2
    local o0
    local o1
    forvalues j = 1/`kx' {
        local o0 `o0' `j'
        local o1 `o1' `=`kx1'+`j''
    }
    local o0 `o0' `=2*`kx1'+`kz1'+1' `kx1'
    local o1 `o1' `=2*`kx1'+`kz1'+2' `=2*`kx1''
    local n0 = `kx' + 2
    forvalues a = 1/`n0' {
        local ia : word `a' of `o0'
        local ja : word `a' of `o1'
        matrix `B'[1,`ia'] = `b0'[1,`a']
        matrix `B'[1,`ja'] = `b1'[1,`a']
        forvalues c = 1/`n0' {
            local ic : word `c' of `o0'
            local jc : word `c' of `o1'
            matrix `VV'[`ia',`ic'] = `V0'[`a',`c']
            matrix `VV'[`ja',`jc'] = `V1'[`a',`c']
        }
    }
    * probit block: Z..., _cons -> positions 2*kx1+1 .. 2*kx1+kz1
    forvalues a = 1/`kz1' {
        local ia = 2*`kx1' + `a'
        matrix `B'[1,`ia'] = `bg'[1,`a']
        forvalues c = 1/`kz1' {
            local ic = 2*`kx1' + `c'
            matrix `VV'[`ia',`ic'] = `Vg'[`a',`c']
        }
    }
    * names
    local names
    foreach v of local xlist {
        local names `names' `eq0':`v'
    }
    local names `names' `eq0':_cons
    foreach v of local xlist {
        local names `names' `eq1':`v'
    }
    local names `names' `eq1':_cons
    foreach v of local zlist {
        local names `names' select:`v'
    }
    local names `names' select:_cons rs0:_cons rs1:_cons
    matrix colnames `B'  = `names'
    matrix colnames `VV' = `names'
    matrix rownames `VV' = `names'
    return matrix b = `B'
    return matrix V = `VV'
end
