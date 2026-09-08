* msat 1.3.0 -- twostep and semipar
* prerequisites: msat.ado 1.3.1, _msat_one.ado 1.2.0, _msat_pwr.ado in the adopath, then . discard
set more off

* ---- Jeu normal (esr_sim) : FIML et deux-etapes doivent coincider ----------
use "esr_sim.dta", clear
movestay y x, select(d = x z)
msat d, mte(.1 .25 .5 .75 .9)
msat d, twostep mte(.1 .25 .5 .75 .9) semipar graph
* attendu twostep : ATT ~ 1.18, ATU ~ -0.34, kappa ~ 1.09 (FIML 1.179 / -0.344 / 1.089)
* attendu semipar : colonne PWR proche de la colonne parametrique a tous les tau

* ---- erreurs asymetriques : FIML casse, deux-etapes non ------------------
use "esr_skew.dta", clear
movestay y x, select(d = x z)
msat d, mte(.1 .25 .5 .75 .9)
* attendu FIML : ATT ~ 2.33, ATU ~ -1.63, kappa ~ 2.82
msat d, twostep mte(.1 .25 .5 .75 .9) semipar graph
* attendu twostep : ATT 1.2588, ATU -0.3449, kappa 1.1669  (= two-step a la main du 8 sept.)
* attendu semipar : PWR proche de la ligne parametrique deux-etapes ; ecart croissant vers les queues avec la ligne FIML
bootstrap att=r(att) atu=r(atu) ate=r(ate) kappa=r(kappa), reps(100) seed(20260908) nodots: ///
    _msat_one y x, select(d = x z) twostep

* ---- gain non lineaire en u : les deux droites perdent la courbe -----------
use "esr_nonlin.dta", clear
movestay y x, select(d = x z)
msat d, twostep mte(.1 .25 .5 .75 .9) semipar graph
* attendu twostep : ATT 1.0880, ATU -0.2101, kappa 0.8782
* attendu semipar : PWR non monotone (0.5 -> 0.7 -> -0.3 environ), la ligne parametrique droite
