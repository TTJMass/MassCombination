# Global impacts on the top mass — implementation plan (Tier 1 + Tier 2)

Status: **planning + hypothesis VERIFIED**, nothing implemented yet. Written 2026-07-23.
Resume cold from this file. Verification artifacts live in `global_impacts_dev/`.

## Context / references
- Reinterpretation top-mass analysis: two-stage design.
  - **Stage 1** = pyconvino cross-section combination (`pyconvino/`), fits combined
    xsec `x` + all ~600 experimental nuisances `alpha`, profiles the nuisances out,
    exports combined values + covariance (+ impacts) as `*_result.npz`.
  - **Stage 2** = mtpole-ttj-pyconvino chi2 mass fit (`mtpole-ttj-pyconvino/fit_object.py`),
    fits the mass `m` (+ PDF nuisances + externalised scale) to `x` using the combined
    covariance. Experimental nuisances are NOT parameters of stage 2.
- Target method: arXiv:2307.04007 "global impacts", as summarised in
  `~/Downloads/CMSWeek-GlobalImpacts.pdf` (A. Gilbert, 13 Feb 2026), slides 4-6.
  - Ordinary/nominal impact (slide 4): fix NP to postfit +-1sigma, refit;
    Gaussian limit = `sigma_mu * rho(mu, nu_i)`.
  - Global impact (slide 5/6): shift the *auxiliary/global observable* by the
    **prefit** sigma, NP floats; Gaussian limit = ordinary * `sigma_post/sigma_prefit`.
  - Key virtue (slide 6, paper Eq 34): global impacts are **additive** —
    `sigma_syst^2 = sum_i (global_impact_i)^2`, valid to sum into groups, because
    aux observables are uncorrelated by construction. Freeze/conditional impacts
    (slide 3) are NOT additive.

## Why two-stage is the right architecture (already established)
The mass is, in the Gaussian/Fisher limit that every impact method assumes, a
linear functional of the combined xsec: `m_hat - m0 = w^T (x - x0)` with
`w = _fisher_weights` (`fit_object.py:260`). For any covariance component `C^(a)`
of `x`, `sigma_a(m)^2 = w^T C^(a) w`. The combined xsec + full covariance are a
sufficient statistic (BLUE/Gauss-Markov); the combination is mass-independent, so
two-stage loses nothing on central value / total unc / stat-syst split. A single
simultaneous fit would only add in-situ constraint of experimental nuisances by
the theory — undesirable (combination is meant theory-agnostic). Keep two-stage.

## Numerical facts VERIFIED (against real NPZs)
Verified on `pyconvino/out/ATLAS813CMS13_corrV2_result.npz` (19 bins, 614 nuis, 48 groups)
and, for cov_full, on a fresh CMSOnly combine (`global_impacts_dev/cmsonly_test.npz`,
348 nuis, 4 bins). Repro scripts noted below.

1. `stat_only_covariance + total_syst_covariance == combined_covariance` to 5e-17.
   => mass stat/syst split is EXACT: `w^T Vstat w + w^T Vsyst w = w^T Ctot w`.
2. **`Vsyst == A · Sigma_post · A^T`** where `A = impact_matrix` (nest x nsys) and
   `Sigma_post` = nuisance-nuisance block of `cov_full`:
   **trace ratio 1.000000, max element-wise err 4.6%** (CMSOnly). The 4.6% is a
   non-Gaussian/asymmetric-response + Hessian-approx artifact, not structural.
3. Diagonal shortcut `A · diag(constraints^2) · A^T`: **67% error** — the exported
   `constraints` (= diag sqrt of Sigma_post, confirmed to 2e-16) is NOT sufficient.
   `Sigma_post` off-diagonal/diagonal ratio = 0.303 (nuisances ~30% correlated postfit).
4. `impact_cov_groups__*` (pyconvino's per-group freeze covariances) do NOT sum to
   Vsyst — overcount ~77x on corrV2 (huge cancelling cross-experiment correlations).
   Confirms they are the freeze/conditional method (slide 3), non-additive.
5. An additive factorization exists and closes: eigen-factor `Sigma_post = L L^T`,
   `V_syst = (A L)(A L)^T = sum_k v_k v_k^T`, and `sum_k (w·v_k)^2 = w^T Vsyst w`
   (0.05% residual). The open piece is mapping additive directions -> physical
   groups, which needs the PREFIT/prior (aux-obs) basis, owned by pyconvino.

### Conclusion enabling the work
The only thing missing to compute genuine additive global impacts on the mass is
the per-source decomposition in the prefit/aux-obs basis. pyconvino has everything
(full Hessian, prior correlation `C_exact`/`pre_sys_corr`). Either export `cov_full`'s
Sigma_post block (cheapest, enables the `A Sigma_post A^T` route) OR — cleaner —
have pyconvino compute per-group additive global-impact VECTORS `v_g` on the
combined observables (paper Eq 30/32) with verified `sum_g v_g v_g^T = Vsyst`.

---

## TIER 1 — exact stat/syst split of the mass (no pyconvino change)
Small, self-contained, uses only arrays already in the NPZ. Do first.

Where: `mtpole-ttj-pyconvino/fit_object.py`, method `doExpStatSystBreakdown`
(currently ~line 1897: profiled refit with stat-only cov + quadrature subtract —
an approximation).

Change: compute the split as an exact linear functional anchored to the existing
profiled `exp_unc` (which correctly captures theory-interp nonlinearity):
```
w  = self._fisher_weights()  -> mass sensitivity vector (already exists)
# In NPZ (normalised fit space) need Vstat, Vsyst mapped through the same
# normalisation Jacobian J used elsewhere (see _build_norm_jacobian / w_abs in
# doUncertaintyBudget: w_abs = J^T @ w). Use w_abs against the ABSOLUTE-space
# stat_only_covariance / total_syst_covariance from the NPZ, consistent with
# how doUncertaintyBudget already contracts.
f          = (w_abs @ Vsyst @ w_abs) / (w_abs @ Ctot @ w_abs)   # exact syst fraction
sigma_syst = self.exp_unc * sqrt(f)
sigma_stat = self.exp_unc * sqrt(1 - f)
```
Guaranteed `sigma_stat^2 + sigma_syst^2 = exp_unc^2`. Add a closure assert.
NOTE: this splits the EXPERIMENTAL uncertainty only; PDF/interp/scale are theory
and already handled separately+correctly (doBreakdown / doScaleVariations).
Verify: run one NPZ fit (e.g. CMS_13TeV_npz stripper), check the printed
stat/syst vs the old method and confirm closure.

---

## TIER 2 — additive global impacts on the mass (needs a pyconvino export)
Two sub-parts.

### 2a. pyconvino: emit additive per-group (and optionally per-nuisance) global-impact vectors
- In `pyconvino/combiner.py` (post-fit section, ~step 12/13, has `H_fit`,
  `impact_matrix`, prior correlation, group defs `config.impact_groups`).
- Implement paper Eq 30/32: for each uncorrelated aux-obs / prior eigen-direction,
  the response of the combined observables `x` when that aux-obs is shifted by its
  prefit sigma and all NPs re-profile. Aggregate per user group -> `v_g` (nest-vector).
  DESIGN DECISION to resolve first: are priors diagonal (independent aux obs ->
  per-nuisance sources are already orthogonal, group vector = quadrature of members)
  or correlated (`C_exact` non-trivial, e.g. corrV2 -> must work in prior eigenbasis
  and cross-experiment-correlated systematics form joint sources)? Inspect the parsed
  prior matrix in the combiner to decide.
- VERIFY inside pyconvino: `sum_g v_g v_g^T == total_syst_covariance` (expect ~<=5%
  matching the A Sigma_post A^T residual; trace should be ~exact).
- Export via `result.py` to_dict/npz as `global_impact_group__<label>` (nest,) each.
  (Keep existing impact_groups / impact_cov_groups untouched — different method.)

### 2b. mtpole-ttj-pyconvino: contract with mass weight, produce the breakdown + forest plot
- New method e.g. `doGlobalImpacts` (or fold into doUncertaintyBudget): read
  `global_impact_group__*`, compute `dm_g = w_abs @ (J-mapped v_g)`, report per-group
  and grouped-into-supergroups; assert `sqrt(sum_g dm_g^2) == sigma_syst(m)` (Tier 1).
- Forest plot mirroring the paper (per-group global impact, sorted), CSV like the
  existing budget_*.csv. Blinding: impacts are shifts/widths (salt-invariant) -> no
  blinding needed, but route through existing `_blind_value` only if any absolute
  mass is printed.
- Deprecate/relabel the misleading current pieces:
  - `direct_impact_MeV = |w_abs @ A[:,i]|` (doUncertaintyBudget) is NOT a global
    impact and does not close — relabel or remove once 2b lands.
  - `global_impact_MeV = sum_i (dm/dtheta_i * pull_i)` is a signed pull-weighted
    coherent SHIFT, a different diagnostic — keep but rename to avoid the "global
    impact" name clash (e.g. `coherent_pull_shift_MeV`).

---

## Open questions / decisions for next session
1. Prior structure (diagonal vs correlated aux obs) in the combination — determines
   the exact 2a construction. Check the parsed prior in combiner for corrV2.
2. Per-group only, or per-nuisance too? (614 nuisances -> per-group=48 is the
   practical paper deliverable; per-nuisance is a big forest plot.)
3. Is the 4.6% `A Sigma_post A^T` residual acceptable, or do we want pyconvino to
   compute v_g from an exact non-Gaussian aux-obs shift (refit) rather than the
   Hessian linearisation? (Trace is exact; likely fine per-group.)
4. Which fits get the new breakdown — all 7 user-facing fits, or just the headline
   full combination?

---

## 2026-07-23 UPDATE — Tier 2a math SOLVED (supersedes the Tier-2 sketch above)
Derived + numerically verified (CMSOnly, fast `compute_impacts=False` 40s combine).

**The genuine per-nuisance global impact on the mass is `dm_i = Cov(m_hat, nu_i)`.**
- `m_hat = w_abs·x`, so `dm_i = w_abs · Sxa[:,i]` where `Sxa = Cov(combined_x, nuisances)`
  = `cov_full[nest-block, nsys-block]`, shape (nest, nsys).
- Identity confirmed: `Sxa == A · Sigma_post` (5.6e-15).
- For UNIT prefit priors (pyconvino uses these — `constraints` defined with prior
  sigma = 1; verify in `_build_prior`), `dm_i = Cov(m,nu_i)` EXACTLY equals the
  slide-5 Gaussian formula `sigma_mu·rho(mu,nu_i)·sigma_post/sigma_pre`. (That's
  literally why slide 5 calls it the "post-fit nuisance-parameter covariance" method.)
- **Additivity is EXACT** (residual 4e-16): `sum_i dm_i^2 = w_abs·(Sxa Sxa^T)·w_abs`,
  ratio 1.000000 for any w. Per-GROUP = quadrature sum `sqrt(sum_{i in g} dm_i^2)`.
- **Global-impacts syst `V_glob = Sxa Sxa^T` != freeze `V_syst`** (12% element-wise,
  trace 0.90 on CMSOnly). They are DIFFERENT decompositions — exactly the slide-6/7
  warning. So Tier 2 gives its OWN stat/syst split, distinct from Tier 1's freeze split.
  Report BOTH, clearly labelled. Do NOT force Tier 2 to close to freeze-V_syst.

**Discarded dead-ends (documented so nobody retries them):**
- per-nuisance `constraint_i · A[:,i]` (diag Sigma_post): 67% error. NO.
- per-group within-block `A_g Sigma_post[g,g] A_g^T`: mass ratios 0.49–1.54, 31.7%
  matrix residual. NO (post-fit cross-group correlations + 119/348 ungrouped nuis).

### STATUS: Tier 2a IMPLEMENTED + VERIFIED (2026-07-23). Uncommitted.
Edits: `pyconvino/combiner.py` (new `CombinationResult.sys_group_labels` field +
built from `cfg.impact_groups` before the return) and `pyconvino/result.py`
(`to_dict` now exports `x_sys_cov` = cov_full[nsys:,:nsys] and `sys_group`).
Verified on a fresh CMSOnly combine: keys present in to_dict AND round-tripped .npz;
closure exact (per-group quadrature == per-nuisance == w·Sxa·Sxa^T·w, ratio 1.000000).
Nice property: `x_sys_cov` only needs cov_full (step 7), so it is available even
with `--no-impacts` — global impacts do NOT require the expensive per-systematic
freeze loop.
**REMAINING for Tier 2:** (1) Sonnet is doing Tier 2b consumer + Tier 1 (contract
corrected via message). (2) The existing fitmatrix combination NPZs predate this
export and LACK x_sys_cov — regenerate at least the nominal combination NPZ (a full
matrix Condor rerun is a separate cost decision, needs user go-ahead). (3) Re-verify
V_glob additivity + magnitude on a real corrV2 (multi-experiment, correlated priors)
combine once one is regenerated — CMSOnly is single-experiment.

**Revised Tier 2a export (TINY — just add to pyconvino result):**
- `x_sys_cov` : (nest, nsys) = `cov_full[nsys:, :nsys]` = Cov(combined_x, nuisance_i),
  absolute units. THE key new array.
- `sys_group` : (nsys,) str — each nuisance's group label from `config.impact_groups`,
  '' / 'ungrouped' if none. (CMSOnly: 119/348 ungrouped — MUST have an ungrouped bucket.)
- `sys_names` already exported (column order).
Implementation: add `sys_group_labels` to `CombinationResult` in combiner.py (build
from `self.config.impact_groups` + `all_sys_names`), export both in `result.to_dict`.
NOTE: standalone single-experiment NPZs used as mtpole INPUT are per-experiment; the
mass fit reads the COMBINATION npz — make sure x_sys_cov is on the same npz the mass
fit consumes (the combination result, e.g. ATLAS813CMS13_corrV2). nest here = n_total
= len(combined_values) = mtpole's n_total.

**Correlated-prior caveat (corrV2):** with non-diagonal prior correlations, the SUM
`sum_i dm_i^2` and per-GROUP quadrature sums stay exactly additive (outer-product
sum), but the per-NUISANCE attribution is basis-dependent (independent sources are
prior eigen-directions). Fine for per-group reporting; note it if doing per-nuisance
forest plots on corrV2. Also: CMSOnly is single-experiment; still should re-verify
`V_glob` additivity + magnitude on a real corrV2 combine (needs cov_full there too,
which today isn't exported — the x_sys_cov export fixes that).

**Tier 2b consumer (given to the Sonnet agent, corrected):** dm = w_abs @ x_sys_cov;
per-group quadrature (incl 'ungrouped'); syst_global = sqrt(sum dm^2); stat_global =
sqrt(w_abs·Ctot·w_abs − syst_global^2); CSV + forest plot + results.json; assert
sqrt(sum_g impact_g^2) == syst_global.

## Repro / artifacts
- `global_impacts_dev/verify_covfull_hypothesis.py` — runs CMSOnly combine (~1420s!)
  and tests the identities. cov_full is NOT exported, so this re-combines in-process.
- `global_impacts_dev/cmsonly_test.npz` — saved arrays (cov_full, impact_matrix,
  Vstat/Vsyst/Ctot, constraints, sys_names) so the identities can be re-checked
  WITHOUT the 1420s re-combine.
- `global_impacts_dev/cmsonly_covtest_result.log` — the verified numbers.
- Combiner in-process API: `from pyconvino.combiner import Combiner` (run from repo
  root, not from inside pyconvino/); `Combiner.from_config(cfg, prefix=..., compute_impacts=True).combine()`;
  result has `.cov_full`, `.impact_matrix`, `.sys_names`, `.constraints`, `.combined_names`;
  `pyconvino.result.to_dict(result)` gives `total_syst_covariance` etc.
- Small setups for fast iteration: `pyconvino/ConvinoSetups/{CMSOnly,ATLAS8Only}/rho_config.txt`.
  (CMSOnly full-impacts combine is still ~1420s because of the per-systematic frozen
  covariance loop; for pure identity checks reuse the saved npz.)
