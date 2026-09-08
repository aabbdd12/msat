* msat_semipar.do -- courbe MTE semi-parametrique (msat, semipar) sur les trois jeux
* -------------------------------------------------------------------------------
* Prerequis (une seule fois) :
*   1. copier _msat_pwr.ado (dossier command/ de msat_1.3_release.zip) dans le
*      MEME repertoire que msat.ado (par ex. PERSONAL : . sysdir list)
*   2. . discard
*   3. verifier : . which _msat_pwr        (doit afficher le chemin, pas "not found")
* Les fichiers esr_sim.dta, esr_skew.dta, esr_nonlin.dta dans le repertoire courant.
* Sortie : trois tableaux (tau, p_tau, MTE parametrique, MTE PWR, e.-t., h, N_eff,
* largeur) et trois graphiques mte_sim.png, mte_skew.png, mte_nonlin.png.
* Copiez-moi la sortie ecran des trois tableaux semipar (c'est ce qui manque a la note B).

set more off
capture which _msat_pwr
if _rc {
    di as err "_msat_pwr.ado introuvable : copiez-le a cote de msat.ado puis tapez discard"
    exit 111
}
local taus ".1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .85 .9"

* ---- 1. design de base (normal) : la courbe PWR doit suivre la droite ------
use "esr_sim.dta", clear
movestay y x, select(d = x z)
msat d, twostep mte(`taus') semipar graph
matrix list r(mte_sp), format(%8.4f)
graph export "mte_sim.png", name(msat_mte) replace width(1400)

* ---- 2. erreurs asymetriques : PWR suit la droite deux-etapes, pas la FIML --
use "esr_skew.dta", clear
movestay y x, select(d = x z)
msat d, mte(`taus')                        // FIML : droite deformee (kappa ~ 2.8)
msat d, twostep mte(`taus') semipar graph  // deux etapes + PWR
matrix list r(mte_sp), format(%8.4f)
graph export "mte_skew.png", name(msat_mte) replace width(1400)

* ---- 3. gain non lineaire en u : seule la PWR suit la courbure --------------
use "esr_nonlin.dta", clear
movestay y x, select(d = x z)
msat d, twostep mte(`taus') semipar graph
matrix list r(mte_sp), format(%8.4f)
graph export "mte_nonlin.png", name(msat_mte) replace width(1400)
