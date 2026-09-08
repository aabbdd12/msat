* =====================================================================
*  VALIDATION de msat 1.2.1 --  delta contre bootstrap-200, courbes par
*  decile contre les vrais effets individuels, MTE contre la verite
*
*  Prerequis : msat.ado et _msat_one.ado (v1.2) dans l'adopath, puis
*              . discard
*  Deux jeux a parametres CONNUS, produits par esr_fiml_reference.py :
*    esr_annex2.dta  effet homogene = 2   (ATT = ATU = ATE = 2.000, kappa = 0)
*    esr_sim.dta     kappa = 1.08         (ATT 1.157 / ATU -0.313 / ATE 0.509)
*
*  Valeurs attendues (Python FIML independant, variance par hessienne ;
*  attendre une concordance a ~10 % sur les SE, pas exacte) :
*
*    esr_sim          estim.     SE delta
*       ATT           1.1792     0.0334
*       ATU          -0.3435     0.0370
*       ATE           0.5081     0.0268
*       kappa         1.0892     0.0373      (verite 1.08)
*       support commun de P(Z) : [0.0213, 0.9910]
*
*    esr_sim, deciles de x  (traites / non traites / tous, puis vrai delta)
*       q1   -0.160  -1.477  -0.915   |  -0.210  -1.429  -0.909
*       q5    0.987  -0.291   0.393   |   0.944  -0.290   0.371
*       q10   2.350   0.957   1.934   |   2.286   1.038   1.913
*
*    esr_sim, deciles de P(Z)
*       q1    2.075  -0.092   0.085   |   1.848  -0.037   0.117
*       q5    1.249  -0.477   0.475   |   1.188  -0.421   0.466
*       q10   1.020  -1.291   0.934   |   1.040  -1.262   0.955
*
*    esr_sim, MTE(u) = 0.5088 + 1.0892*invnormal(1-u)   (verite 0.506 + 1.08*.)
*       u=.05  2.300 (0.067)   u=.25  1.244 (0.037)   u=.50  0.509 (0.027)
*       u=.75 -0.226 (0.037)   u=.95 -1.283 (0.067)
* =====================================================================
set more off
clear all
discard

* ---------------------------------------------------------------------
* JEU 1 : esr_annex2.dta  --  effet homogene, kappa = 0
* ---------------------------------------------------------------------
use "esr_annex2.dta", clear
summarize delta_true                       // 2.000 exactement

xi: movestay income educ , select(treatment i.region inst)

* (a) methode delta -- defaut ; kappa doit etre non significatif
msat treatment
di as txt _n "  -> attendu : ATT ~ 1.981 (SE ~0.030), kappa ~ 0 (non significatif)"

* (b) estimations ponctuelles seules (chemin rapide utilise par le bootstrap)
msat treatment, nose

* (c) bootstrap-200, toutes les etapes refaites a chaque tirage
bootstrap att=r(att) atu=r(atu) ate=r(ate) kappa=r(kappa), reps(200) ///
    seed(20260907) nodots: ///
    _msat_one income educ, select(treatment i.region inst) xi
di as txt _n "  -> le SE bootstrap doit etre du meme ordre que le SE delta (~0.03)"


* ---------------------------------------------------------------------
* JEU 2 : esr_sim.dta  --  kappa = 1.08, ATT != ATU
* ---------------------------------------------------------------------
use "esr_sim.dta", clear
summarize delta_true if d==1               // ATT vrai  1.157
summarize delta_true if d==0               // ATU vrai -0.313
summarize delta_true                       // ATE vrai  0.509

movestay y x, select(d = x z)

* (a) effets, kappa, courbe par decile de x, MTE, graphiques
msat d, generate(dhat) pscore(pz) rank(x) nq(10) ///
        mte(.05 .1 .25 .5 .75 .9 .95) graph
di as txt _n "  -> attendu : ATT ~ 1.179 (0.033), ATU ~ -0.344 (0.037), kappa ~ 1.089 (0.037)"

* verification directe : effet attendu vs effet vrai par decile de x
xtile qx = x, nq(10)
tabulate qx d, summarize(dhat) means nostandard nofreq
tabulate qx d, summarize(delta_true) means nostandard nofreq

* (b) courbe par decile du score de participation : ATT(a) et ATU(a)
*     doivent se deformer en sens opposes (kappa > 0)
msat d, rank(pz) nq(10) graph
xtile qp = pz, nq(10)
tabulate qp d, summarize(dhat) means nostandard nofreq
tabulate qp d, summarize(delta_true) means nostandard nofreq

* (c) point estimates only + bootstrap-200 (kappa inclus)
msat d, nose
bootstrap att=r(att) atu=r(atu) ate=r(ate) kappa=r(kappa), reps(200) ///
    seed(20260907) nodots: ///
    _msat_one y x, select(d = x z)


* ---------------------------------------------------------------------
* Lecture
* ---------------------------------------------------------------------
* - Estimations identiques a 1e-4 pres entre msat et Python : formules OK.
* - SE delta ~ SE bootstrap ~ valeurs attendues : variance OK.
* - Par decile, dhat moyen et delta_true moyen doivent differer de moins de
*   ~2 SE (delta_true contient en plus la composante orthogonale du gain,
*   d'ecart-type ~1.55, moyennee sur ~2000 obs par decile : ~0.035).
* - MTE : droite de pente kappa en invnormal(1-u) ; drapeau "observed" sur
*   [0.02, 0.99] environ.
