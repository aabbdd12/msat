# msat — treatment effects after endogenous switching regression

`msat` is a Stata post-estimation command for `movestay` (Lokshin and Sajaia 2004).
After an endogenous switching regression it returns

- the four conditional expectations E[Y_j | D] and the ATT, ATU and ATE with delta-method standard errors (parameter uncertainty plus the sampling component);
- the selection-on-gains parameter **kappa = rho1·sigma1 − rho0·sigma0**, with its standard error and Wald test — ATT − ATU = kappa·(mean lambda1 + mean lambda0);
- the expected-effect curves by quantile of any ranking variable (`rank()`, `nq()`, `graph`);
- the marginal treatment effect at chosen percentiles of the participation unobservable, with common-support flags (`mte()`);
- a **two-step** re-estimation of the model — probit, then OLS on X and the Mills ratio in each regime — which needs only a normal participation error and a linear conditional mean, and stays exact when the likelihood breaks under skewed regime errors (`twostep`);
- a **semiparametric MTE** by percentile-weights regression on the probit score, with a derivative-aware bandwidth, next to the parametric line (`semipar`).

Survey weights (`svyset`) and household size (`hsize()`) enter every average.

## Installation

From Stata 12 or later:

```stata
ssc install movestay                     // required
net install msat, from("https://raw.githubusercontent.com/aabbdd12/msat/main/") replace
net get msat, from("https://raw.githubusercontent.com/aabbdd12/msat/main/")   // validation datasets and do-files (optional)
discard
help msat
```

## Example

```stata
use esr_sim, clear                       // fetched by net get (see Installation)
movestay y x, select(d = x z)
msat d                                   // ATT, ATU, ATE, kappa
msat d, rank(x) nq(10) graph             // effect curves by decile of x
msat d, mte(.1 .25 .5 .75 .9)            // MTE with support flags
msat d, twostep mte(.1(.05).9) semipar graph   // two-step route + semiparametric MTE
```

Exact standard errors for the two-step route, by bootstrap over the whole procedure:

```stata
bootstrap att=r(att) atu=r(atu) ate=r(ate) kappa=r(kappa), reps(200): ///
    _msat_one y x, select(d = x z) twostep
```

## Files

| File | Role |
|---|---|
| `msat.ado`, `msat.sthlp` | the command and its help file |
| `_msat_one.ado` | bootstrap helper: re-estimates `movestay` and recomputes the effects on each resample |
| `_msat_pwr.ado` | private percentile-weights-regression engine used by `semipar` (derived from `gepwreg` 1.4; not meant to be called directly) |
| `esr_*.dta` | four simulated datasets with known parameters (see below) |
| `msat_validate.do`, `msat13_test.do`, `msat_semipar.do` | do-files that reproduce the validation checks on those datasets |
| `msat_note.pdf` | the technical note: model, estimands, the two routes, standard errors, heterogeneity, Monte Carlo and illustrations |

## Validation datasets

| Dataset | n | Design | True ATT / ATU / kappa |
|---|---|---|---|
| `esr_sim.dta` | 20,000 | joint normality holds, kappa = 1.08 | 1.157 / −0.313 / 1.08 (sample) |
| `esr_annex2.dta` | 5,000 | homogeneous effect, kappa = 0 (variables `income treatment educ inst region`) | 2 / 2 / 0 |
| `esr_skew.dta` | 10,000 | skewed regime errors: FIML breaks, two-step does not | 1.148 / −0.329 / 1.08 |
| `esr_nonlin.dta` | 10,000 | gain nonlinear in u: both lines lose the curve, `semipar` follows it | 1.033 / −0.174 / 1.08 (projection) |

## Citation

Araar, A. (2026). *Treatment Effects after Endogenous Switching Regression: the msat Command.* PEP technical note (msat_note.pdf in this repository). Software: https://github.com/aabbdd12/msat.

The semiparametric engine implements the percentile-weights regression of Araar, A. (2026), *Exploring Heterogeneous Effects: Quantile Models and Percentile Weights Regression*, Zenodo, https://doi.org/10.5281/zenodo.20315684 (Stata command `gepwreg`).

## References

Lokshin, M., and Z. Sajaia (2004). Maximum likelihood estimation of endogenous switching regression models. *Stata Journal* 4(3): 282–289.

Heckman, J. J., and E. Vytlacil (2005). Structural equations, treatment effects, and econometric policy evaluation. *Econometrica* 73(3): 669–738.

## License

MIT — see `LICENSE`.
