{smcl}
{* *! version 1.3.2  09sep2026}{...}
{vieweralsosee "movestay" "help movestay"}{...}
{vieweralsosee "mspredict" "help mspredict"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "etregress" "help etregress"}{...}
{vieweralsosee "teffects" "help teffects"}{...}
{vieweralsosee "bootstrap" "help bootstrap"}{...}
{viewerjumpto "Syntax" "msat##syntax"}{...}
{viewerjumpto "Description" "msat##description"}{...}
{viewerjumpto "Options" "msat##options"}{...}
{viewerjumpto "Remarks" "msat##remarks"}{...}
{viewerjumpto "Bootstrap" "msat##bootstrap"}{...}
{viewerjumpto "Examples" "msat##examples"}{...}
{viewerjumpto "Stored results" "msat##results"}{...}
{viewerjumpto "References" "msat##references"}{...}
{viewerjumpto "Author" "msat##author"}{...}
{title:Title}

{phang}
{bf:msat} {hline 2} Treatment effects after endogenous switching regression:
ATT, ATU, ATE, selection on gains (kappa), effect curves by quantile of a
ranking variable, and marginal treatment effects (parametric and
semiparametric), by full-information or two-step estimation, with
delta-method standard errors


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:msat} {it:treatvar} {ifin} [{cmd:,} {it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt hs:ize(varname)}}weight each observation by {it:varname}
(e.g. household size){p_end}
{synopt:{opt exp:and(yes)}}outcome was modelled in logs; report effects on the
level of the outcome{p_end}

{syntab:Estimation}
{synopt:{opt two:step}}re-estimate the switching regression in two steps
(probit, then OLS on X and the Mills ratio in each regime) from
{cmd:movestay}'s variable lists, and report the effects from these
estimates; robust to non-normal outcome errors{p_end}
{synopt:{opt dep:var(varname)}}outcome variable, if {cmd:msat} cannot recover
it from {cmd:movestay}'s equation names{p_end}
{synopt:{opt semi:par}}with {opt mte()}: add the semiparametric MTE by
percentile-weights regression on the probit score (engine {cmd:_msat_pwr},
shipped with {cmd:msat}){p_end}

{syntab:SE}
{synopt:{opt nose}}report point estimates only; skips the delta-method
computation (for use inside {helpb bootstrap}){p_end}
{synopt:{opt step(#)}}relative step of the numerical gradient; default is
{cmd:step(0.0001)}{p_end}

{syntab:Heterogeneity}
{synopt:{opt gen:erate(newvar)}}save the individual expected effects
E[Y1-Y0 | D, X, Z]{p_end}
{synopt:{opt rank(varname)}}report the expected effect by quantile of
{it:varname}, for the treated, the untreated and all{p_end}
{synopt:{opt nq(#)}}number of quantile groups for {opt rank()}; default is
{cmd:nq(10)}{p_end}
{synopt:{opt mte(numlist)}}report the marginal treatment effect at the listed
percentiles (strictly between 0 and 1) of the participation unobservable{p_end}
{synopt:{opt ps:core(newvar)}}save the estimated probability of treatment
P(D=1 | Z){p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt gr:aph}}draw the {opt rank()} curve and/or the {opt mte()} curve
with confidence bands{p_end}
{synopt:{opt gropt:s(twoway_options)}}options passed to {helpb twoway}{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
{cmd:msat} is a postestimation command for {helpb movestay}. {it:treatvar}
is the binary switching variable used in {cmd:movestay}'s {cmd:select()}
option. The estimation sample is {cmd:e(sample)} intersected with the
{it:if} and {it:in} qualifiers.


{marker description}{...}
{title:Description}

{pstd}
{cmd:msat} computes the average treatment effect on the treated (ATT), on the
untreated (ATU), and the population average treatment effect (ATE) from the
parameters of an endogenous switching regression fitted by {cmd:movestay}
(Lokshin and Sajaia 2004). The effects are averages of individual differences
between the conditional expectations of the two regimes, each evaluated on the
subsample where the individual is actually observed.

{pstd}
Standard errors account for two sources of uncertainty: the sampling variability
of the averaged individual effects (design-based, honouring {cmd:svyset}), and
the estimation uncertainty of the {cmd:movestay} parameters, propagated by the
delta method through {cmd:e(V)}.

{pstd}
The parameter kappa = rho1*sigma1 - rho0*sigma0, the covariance between the
selection error and the unobserved individual gain, is reported with its
standard error; its z-statistic is a Wald test of selection on gains, which is
what distinguishes the switching regression from the single-equation
endogenous treatment model ({helpb etregress}), where kappa = 0 by
construction.

{pstd}
Two further outputs describe the heterogeneity of the effect. {opt rank()}
averages the individual expected effects within quantile groups of any
exogenous ranking variable (an initial endowment, a welfare measure, or the
propensity score itself), separately for the treated and the untreated; when
the ranking variable is the selection index, the two curves deform in opposite
directions if and only if kappa is not zero. {opt mte()} evaluates the
marginal treatment effect MTE(u) = E[X]'(b1-b0) + kappa*invnormal(1-u) at
chosen percentiles u of the participation unobservable (u close to 0: the
individuals most eager to participate), and flags whether each u lies inside
the common support of the estimated probability of treatment or is an
extrapolation of the normal model.

{pstd}
{cmd:msat} is the single entry point after {cmd:movestay}. Its default path
uses {cmd:movestay}'s full-information estimates. With {opt twostep} it
re-estimates the same model by probit followed by least squares of each regime
on X and its Mills ratio, using the weights and variable lists of the
{cmd:movestay} call, and reports the same table from these estimates. The
two-step estimator needs only a normal participation error and a conditional
mean of the regime errors linear in it, where full-information maximum
likelihood needs the entire joint law: on designs with skewed outcome errors
the likelihood fails (correlations driven to +/-1, kappa several times its
value) and the two-step estimator does not. Running both and comparing them
is the practical test of joint normality. Under {opt twostep} the propensity
score, the common support and the MTE flags come from the probit.

{pstd}
With {opt mte()} and {opt semipar}, {cmd:msat} runs a percentile-weights
regression (Araar 2016, 2023; private engine {cmd:_msat_pwr}, derived from
{cmd:gepwreg} and independent of the user's installation of it) ranked on the
probit score, with a derivative-aware bandwidth, at each listed percentile,
and prints the
semiparametric MTE next to the parametric line evaluated at the same score
value. A gap between the two that grows towards the tails is a departure
from joint normality made visible; the width in p of each kernel window is
the support statistic of the semiparametric estimate.

{pstd}
All weighted averages use the survey weight declared through {cmd:svyset}, if
any, multiplied by {opt hsize()}.

{pstd}
{cmd:msat} restores {cmd:movestay}'s estimation results on exit, so it can be
called repeatedly and other postestimation commands remain available.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt hsize(varname)} weights each observation by {it:varname} when averaging
the individual effects, in addition to any survey weights declared through
{cmd:svyset}. Typical use: household size, to obtain individual-level rather
than household-level averages.

{phang}
{opt expand(yes)} specifies that {cmd:movestay} was fitted on the logarithm of
the outcome and that effects are wanted on the outcome itself. Each conditional
expectation is exponentiated before differencing. Note that exp(E[ln Y]) is the
conditional geometric mean, not E[Y]; for a retransformation to the arithmetic
mean see the smearing estimator of Duan (1983). Standard errors are propagated
through the exponential by the delta method.

{dlgtab:Estimation}

{phang}
{opt twostep} replaces {cmd:movestay}'s estimates by the two-step ones. The
outcome and the regressors of the regime equations, the selection regressors
and the weights are read from {cmd:movestay}'s results; the regime
regressions use robust variances. The parameter vector is laid out as
{cmd:movestay}'s, with {cmd:rs0:_cons} and {cmd:rs1:_cons} (the coefficients
on the Mills ratio, equal to rho_j*sigma_j) in place of the variance and
correlation parameters, so that kappa = rs1 - rs0. The delta-method variance
treats the three estimations as independent (block-diagonal), which ignores
the covariance induced by the generated regressor; use
{cmd:bootstrap ... : _msat_one ..., twostep} for exact standard errors.

{phang}
{opt depvar(varname)} names the outcome when it cannot be recovered
({cmd:movestay} names the regime equations {it:depvar}0 and {it:depvar}1).

{phang}
{opt semipar} requires {opt mte()}. For each listed percentile u, the
percentile-weights regression of the outcome on the centred covariates, their
interactions with the score, the score P and P^2 is run, ranked on the probit
score at tau = u, with a bandwidth chosen for the derivative in P: pilot
h0 = 1.6*sigma_p*n^(-1/7), curvature from a pilot curve, then the empirical
MSE of the target minimised over a grid of multiples of h0 (capped at 0.30 in
rank space). The MTE is {cmd:_b[P] + 2*q_tau*_b[P^2]} at the tau-quantile
q_tau of the score, with the influence-function standard error of Deville
(1999), which excludes the estimation of the score. Results are stored in
{cmd:r(mte_sp)} and overlaid on the MTE graph. The engine {cmd:_msat_pwr.ado}
ships with {cmd:msat}; the same bandwidth rule is available to
{cmd:gepwreg} users as its {cmd:target()} option (version 1.4).

{dlgtab:SE}

{phang}
By default the variance of each effect is computed as G*e(V)*G' + S, where G
is the numerical gradient of the effect with respect to the full parameter
vector {cmd:e(b)} of {cmd:movestay} (in its estimated parameterization:
ln sigma and atanh rho), e(V) is {cmd:movestay}'s variance matrix, and S is
the design-based sampling variance of the average of the individual effects.
The two components are reported separately below the table and are stored in
{cmd:r()}. The parameter component G*e(V)*G' is usually the dominant one.

{phang}
{opt nose} suppresses the standard errors and reports point estimates only.
The delta-method gradient is not computed and {cmd:svy: ratio} is not called,
which makes the command about twice as fast; this is the path used by
{cmd:_msat_one} inside {helpb bootstrap}
(see {help msat##bootstrap:Bootstrap} below).

{phang}
{opt step(#)} sets the relative step of the central-difference gradient,
h_j = # * max(1, |b_j|). The default 0.0001 is adequate for all cases tested.
Increase it if {cmd:movestay} converged with a very flat likelihood.

{dlgtab:Heterogeneity}

{phang}
{opt generate(newvar)} saves, for each observation of the estimation sample,
the expected effect d_i = E[Y1 - Y0 | D_i, X_i, Z_i] whose weighted averages
are the ATT, ATU and ATE. Use it with any external tool to describe the
distribution of expected gains; do not rank it on the outcome or on a
consequence of the treatment.

{phang}
{opt rank(varname)} reports the weighted mean of d_i within each quantile
group of {it:varname} (groups formed on the whole estimation sample), for the
treated, the untreated, and all observations, with standard errors that
combine the delta-method and sampling components exactly as for the aggregate
effects. Empty cells are reported as missing. Rank on a pre-treatment or
exogenous variable.

{phang}
{opt nq(#)} sets the number of quantile groups; default 10.

{phang}
{opt pscore(newvar)} saves the estimated probability of treatment Phi(Z'g),
which can then be used as the ranking variable of a second call.

{phang}
{opt mte(numlist)} reports MTE(u) at each listed percentile u in (0,1). Under
the joint normality of {cmd:movestay}, MTE(u) is linear in invnormal(1-u) with
slope kappa; its standard error is obtained from the joint delta-method
variance of E[X]'(b1-b0) and kappa. The support flag compares u with the
common support of the estimated probability of treatment, [max of the two
group minima, min of the two group maxima]. With {opt expand(yes)} the MTE is
reported on the scale of the linear index.

{dlgtab:Reporting}

{phang}
{opt level(#)}; see {helpb estimation options##level():[R] Estimation options}.
Confidence intervals are normal-based.

{phang}
{opt graph} draws the {opt rank()} curve (graph name {cmd:msat_rank}) and/or
the MTE curve on the grid u = 0.01(0.01)0.99 with the common support marked by
vertical lines (graph name {cmd:msat_mte}).

{phang}
{opt gropts(twoway_options)} are appended to the {cmd:twoway} call.


{marker remarks}{...}
{title:Remarks}

{pstd}
{ul:Conditional expectations and the sign of the Mills ratio}

{pstd}
With the switching equation D = 1[Z'g + u > 0], u ~ N(0,1), and regime
equations Y_j = X'b_j + e_j, corr(e_j, u) = rho_j, sd(e_j) = sigma_j,
the four conditional expectations are (Lokshin and Sajaia 2004, eqs 5-8):

{p 8 12 2}E[Y1 | D=1] = X'b1 + rho1*sigma1*lambda1{p_end}
{p 8 12 2}E[Y0 | D=1] = X'b0 + rho0*sigma0*lambda1{p_end}
{p 8 12 2}E[Y1 | D=0] = X'b1 - rho1*sigma1*lambda0{p_end}
{p 8 12 2}E[Y0 | D=0] = X'b0 - rho0*sigma0*lambda0{p_end}

{pstd}
where lambda1 = phi(Z'g)/Phi(Z'g) and lambda0 = phi(Z'g)/(1-Phi(Z'g)).
The sign of the correction is fixed by the conditioning subsample, +lambda1
on the treated and -lambda0 on the untreated, whichever regime is being
predicted.

{pstd}
The effects are

{p 8 12 2}ATT = mean over D=1 of [ X'(b1-b0) + kappa*lambda1 ]{p_end}
{p 8 12 2}ATU = mean over D=0 of [ X'(b1-b0) - kappa*lambda0 ]{p_end}
{p 8 12 2}ATE = p*ATT + (1-p)*ATU = mean over all of X'(b1-b0){p_end}

{pstd}
with kappa = rho1*sigma1 - rho0*sigma0. The parameter kappa is the covariance
between the selection error and the unobserved individual gain; its sign
determines the ordering of ATT and ATU (ATT > ATU if and only if kappa > 0),
and the Mills terms cancel exactly in the ATE.

{pstd}
{ul:Standard errors}

{pstd}
The individual effects are smooth functions of {cmd:e(b)}; the delta method
is therefore valid and the bootstrap over the whole procedure (re-estimating
{cmd:movestay} on each resample) is a valid check. On datasets with known
parameters, the delta-method standard errors of version 1.1 agree with a
200-replication bootstrap and with the Monte Carlo standard deviation of the
estimator to within about 8 percent, with 95-percent coverage between 94.5 and
97 percent.

{pstd}
{ul:Common support and heterogeneity}

{pstd}
The counterfactual expectations extrapolate the parametric model into regions
of the selection index where one regime may not be observed. It is good
practice to report the range of the estimated probability of treatment on
which both regimes are populated, and to restrict the estimation sample to it
before running {cmd:movestay}; the ATU is the estimand most exposed to this
extrapolation. {cmd:msat} reports the common support of the estimated
probability of treatment, and the {opt mte()} table says which percentiles of
the participation unobservable are observed and which are extrapolated.

{pstd}
Three kinds of heterogeneity should be kept apart: the observable one,
X'(b1-b0), which varies with X; the unobservable one correlated with
participation, kappa*lambda, which is what {opt rank()} on the selection index
and {opt mte()} display; and the unobservable one orthogonal to participation,
whose variance is not identified, so that quantiles of the individual gain
Y1 - Y0 are not identified without further assumptions.


{marker bootstrap}{...}
{title:Bootstrap}

{pstd}
The companion program {cmd:_msat_one} re-estimates {cmd:movestay} and
recomputes the effects on each resample, so that every data-dependent step is
inside the bootstrap loop. It mirrors {cmd:movestay}'s syntax:

{p 8 17 2}
{cmd:bootstrap} {cmd:att=r(att) atu=r(atu) ate=r(ate)}{cmd:,} {opt reps(#)}
[{it:bootstrap_options}]{cmd::} {cmd:_msat_one} {it:depvar} {it:indepvars}
{ifin}{cmd:,} {opt sel:ect(treatvar [=] varlist_s)} [{opt xi}
{opt hs:ize(varname)} {opt exp:and(yes)}]

{pstd}
Specify {opt xi} if {cmd:movestay} is to be run with the {cmd:xi:} prefix.
With {cmd:cluster()} or {cmd:strata()} in the {cmd:bootstrap} options, the
resampling honours the survey design. The bootstrap is slower than the delta
method by a factor equal to the number of replications and is intended as a
validation device, not as the routine variance estimator.


{marker examples}{...}
{title:Examples}

{pstd}Homogeneous-effect validation data (true ATT = ATU = ATE = 2){p_end}
{phang2}{cmd:. use esr_annex2, clear}{p_end}
{phang2}{cmd:. xi: movestay income educ, select(treatment i.region inst)}{p_end}
{phang2}{cmd:. msat treatment}{p_end}

{pstd}Point estimates only{p_end}
{phang2}{cmd:. msat treatment, nose}{p_end}

{pstd}Bootstrap check with full re-estimation{p_end}
{phang2}{cmd:. bootstrap att=r(att) atu=r(atu) ate=r(ate), reps(200) seed(1): ///}{p_end}
{phang2}{cmd:      _msat_one income educ, select(treatment i.region inst) xi}{p_end}

{pstd}Outcome in logs, effects on levels, household-size weighting{p_end}
{phang2}{cmd:. movestay lnexp educ age, select(program = educ age dist)}{p_end}
{phang2}{cmd:. msat program, expand(yes) hsize(hhsize)}{p_end}

{pstd}Selection on gains (true kappa = 1.08): kappa with its Wald test,
effect curves by decile of x and of the propensity score, MTE at chosen
percentiles, both graphs{p_end}
{phang2}{cmd:. use esr_sim, clear}{p_end}
{phang2}{cmd:. movestay y x, select(d = x z)}{p_end}
{phang2}{cmd:. msat d, generate(dhat) pscore(pz) rank(x) nq(10) mte(.05 .1 .25 .5 .75 .9 .95) graph}{p_end}
{phang2}{cmd:. msat d, rank(pz) graph}{p_end}

{pstd}Two-step estimates and the semiparametric MTE side by side{p_end}
{phang2}{cmd:. movestay y x, select(d = x z)}{p_end}
{phang2}{cmd:. msat d, mte(.1 .25 .5 .75 .9) semipar graph}{p_end}
{phang2}{cmd:. msat d, twostep mte(.1 .25 .5 .75 .9) semipar graph}{p_end}
{phang2}{cmd:. bootstrap att=r(att) atu=r(atu) ate=r(ate) kappa=r(kappa), reps(200): ///}{p_end}
{phang2}{cmd:      _msat_one y x, select(d = x z) twostep}{p_end}

{pstd}Survey data: declare the weights before {cmd:movestay}{p_end}
{phang2}{cmd:. svyset psu [pw=wgt], strata(strata)}{p_end}
{phang2}{cmd:. movestay lnexp educ age [pw=wgt], select(program = educ age dist)}{p_end}
{phang2}{cmd:. msat program, expand(yes) hsize(hhsize) rank(land) nq(5)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:msat} stores the following in {cmd:r()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:r(att)}}average treatment effect on the treated{p_end}
{synopt:{cmd:r(atu)}}average treatment effect on the untreated{p_end}
{synopt:{cmd:r(ate)}}average treatment effect{p_end}
{synopt:{cmd:r(se_att)}}standard error of {cmd:r(att)}{p_end}
{synopt:{cmd:r(se_atu)}}standard error of {cmd:r(atu)}{p_end}
{synopt:{cmd:r(se_ate)}}standard error of {cmd:r(ate)}{p_end}
{synopt:{cmd:r(kappa)}}rho1*sigma1 - rho0*sigma0{p_end}
{synopt:{cmd:r(se_kappa)}}standard error of {cmd:r(kappa)}{p_end}
{synopt:{cmd:r(m)}}mean observable gain E[X]'(b1-b0) (equals the ATE
without {opt expand(yes)}){p_end}
{synopt:{cmd:r(se_m)}}standard error of {cmd:r(m)}{p_end}
{synopt:{cmd:r(p_lo)}}, {cmd:r(p_hi)} bounds of the common support of the
estimated probability of treatment{p_end}
{synopt:{cmd:r(vd_att)}}delta-method (parameter) variance component, ATT{p_end}
{synopt:{cmd:r(vd_atu)}}same, ATU{p_end}
{synopt:{cmd:r(vd_ate)}}same, ATE{p_end}
{synopt:{cmd:r(vs_att)}}sampling variance component, ATT{p_end}
{synopt:{cmd:r(vs_atu)}}same, ATU{p_end}
{synopt:{cmd:r(vs_ate)}}same, ATE{p_end}

{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:r(G)}}K x p gradient of all reported quantities with respect to
{cmd:e(b)}; rows 1-5 are ATT, ATU, ATE, m, kappa; not stored with {opt nose}{p_end}
{synopt:{cmd:r(curve)}}with {opt rank()}: nq x 7 matrix (q, att, se_att, atu,
se_atu, all, se_all){p_end}
{synopt:{cmd:r(mte)}}with {opt mte()}: one row per percentile (u, mte, se,
support flag){p_end}
{synopt:{cmd:r(psupport)}}1 x 2 common support of the estimated probability of
treatment{p_end}
{synopt:{cmd:r(mte_sp)}}with {opt semipar}: one row per percentile (tau, p_tau,
parametric MTE at p_tau, PWR estimate, se, h, N_eff, width in p){p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Araar, A. 2015. The treatment effect: comparing the ESR and PSM methods with
an artificial example. PEP technical note.

{phang}
Duan, N. 1983. Smearing estimate: a nonparametric retransformation method.
{it:Journal of the American Statistical Association} 78: 605-610.

{phang}
Heckman, J. J., S. Urzua, and E. Vytlacil. 2006. Understanding instrumental
variables in models with essential heterogeneity. {it:Review of Economics and
Statistics} 88: 389-432.

{phang}
Lokshin, M., and Z. Sajaia. 2004. Maximum likelihood estimation of endogenous
switching regression models. {it:Stata Journal} 4: 282-289.

{phang}
Maddala, G. S. 1983. {it:Limited-Dependent and Qualitative Variables in
Econometrics}. Cambridge: Cambridge University Press.


{marker author}{...}
{title:Author}

{pstd}
Abdelkrim Araar{break}
Universit{c e'} Laval and Partnership for Economic Policy (PEP){break}
abdelkrimaraar@gmail.com

{pstd}
Version 1.3.2, September 2026. {cmd:msat} reads the results of either build
of {cmd:movestay}: the Stata Journal build (st0071_2, version 2.0.0, ancillary
parameters lns1 lns2 r1 r2 with index 1 the treated regime, selection
equation named after the treatment variable) and the later build (lns0 lns1
r0 r1, equation "select"); regime 0 is always the untreated regime in the
output. Results differ from version 1.0 (April 2015):
the two counterfactual expectations now follow Lokshin and Sajaia (2004,
eqs 6 and 7), and standard errors include parameter-estimation uncertainty.
