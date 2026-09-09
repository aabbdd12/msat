# Changelog

## 1.3.2 — 2026-09-09
- Reads both `movestay` builds: Stata Journal st0071_2 (v2.0.0; equations `<y>_1`/`<y>_0`, ancillary `lns1 lns2 r1 r2` with index 1 = treated, selection equation named after the treatment variable) and the later build (`<y>0`/`<y>1`, `lns0 lns1 r0 r1`, equation `select`). Fixes `lns0:_cons not found` after the SJ build.

## 1.3.1 — 2026-09-08
- `semipar`: results kept for the graph (`return matrix mte_sp` moved after the plot; fixes r(111)).
- Private engine renamed `_msat_pwr.ado` (prefix `_mpwr_`), no dependency on the user's `gepwreg`.

## 1.3.0 — 2026-09-08
- `twostep`: probit then OLS on X and the Mills ratio in each regime; same four expectations, kappa and MTE; block-diagonal delta-method SEs (bootstrap with `_msat_one, twostep` for exact SEs).
- `semipar` (with `mte()`): semiparametric MTE by percentile-weights regression on the probit score, derivative-aware bandwidth (n^(-1/7) pilot, empirical MSE search), table and overlay graph, `r(mte_sp)`.
- `depvar()` option.

## 1.2.1 — 2026-09-07
- MTE standard error fixed for u > 0.5 (macro parsing of a negative quantile).

## 1.2.0 — 2026-09-07
- `kappa` with Wald test; `generate()`, `pscore()`, `rank()`/`nq()` expected-effect curves with graph; `mte()` with common-support flags; survey weights (`svyset`) and `hsize()` in every average.

## 1.1.2 — 2026-09
- Delta-method standard errors on the full parameter vector plus the sampling component; validated against a 200-replication bootstrap.
