# `--tnp` / `--expScale` mass-fit extension — implementation plan

Status: **plan approved by user, nothing implemented yet**. Written 2026-08-04.
Resume cold from this file. Full context also in the session's Claude Code plan-mode
file (`~/.claude/plans/giggly-wondering-lake.md`, local to the assistant, not synced —
this file is the durable copy).

## Context

`ttbarj-nnlo-cms-atlas-analysis/convertData.py` (see that repo's own history) now
exports two NNLO STRIPPER-only fields into the theory JSON, computed only for PDF
family `PDF4LHC21_40`, orders LO/NLO/NNLO (not `NNLO_noVVF`):
- `tnp_eigenvectors[family][obs_label]`: theory nuisance parameters (TNPs), a
  Chebychev0-based shape-uncertainty model (`ttbarj-nnlo-cms-atlas-analysis/lib/tnps.py`),
  computed at the fully central "CC" scale point only. Count depends on order: 3/6/9
  eigenvectors for LO/NLO/NNLO. Two-sided and symmetric by construction (the model's
  `get_prediction_next_term` is linear in the eigenvector parameter, so +delta/-delta
  are exact negations of each other).
- `normalized_expanded[family][obs_label]`: the order-by-order Taylor-expanded ratio
  rho = dsigma/sigma (`plot-pQCD-normexp.py`'s definition), at all 7 (muR,muF) points —
  same grid as the classic scale variation, just a different definition of the ratio.

Goal: extend `mtpole-ttj-pyconvino`'s mass fit (`doFit.py` -> `fit_object.py`) with:
- `--tnp`: add TNPs as free, unit-Gaussian-constrained nuisance parameters in the
  Minuit fit (profile-likelihood style, same spirit as PDF eigenvectors, but kept in
  their own separate parameter block -- **not merged into `self.PDF_vars`**, so PDF
  and TNP stay cleanly distinguishable everywhere). Replaces the externalized scale
  envelope entirely (`doScaleVariations` is skipped when `--tnp`).
- `--expScale`: swap the classic same-order (muR,muF) ratio for the expanded-ratio
  definition as the source of the *existing* externalized-envelope machinery
  (`doScaleVariations`/`getScaleVariationsMeas`/`setGlobalRescalingFactorTheory`).
  Still externalized, not a fit nuisance -- just different underlying numbers.
- Both fully correlated across all measurements and center-of-mass energies (one
  shared nuisance parameter / one shared externalized shift, not per-measurement).
- `--tnp` gets a new post-fit constraints/pulls plot, analogous to the existing
  `PDF_postfit` plot in `plotFitResults()`.

## Key existing mechanics (mtpole-ttj-pyconvino/fit_object.py)

- `self.PDF_vars` (shape `(n_pdf, n_bins)`): fractional per-bin shift per PDF
  nuisance, concatenated across `self.exp_input.input_meas` -- one shared free
  Minuit parameter multiplies the same row across all measurements' bins, which is
  what makes a nuisance "fully correlated across measurements."
- **The existing `--fitScales` code path is structurally the template to mirror for
  `--tnp`**: it already demonstrates "PDF block + one more independent nuisance
  block after it" throughout the file:
  - `initialiseMinuit()`: `params = append(mass, zeros(len(PDF_vars)), [0,0] if fitScales)`.
  - `_chi2()`: `priors = sum(rescaled_PDF_priors**2)`, `+= params[-2]**2 + params[-1]**2 if fitScales`.
  - `getInterpolatedXsec()`'s non-`_mdv` branch: `PDF_weights = ...` then a separate
    `if self.fitScales: scale_weights = array([mur_var*params[-2], muf_var*params[-1]])`,
    multiplied in as its own factor.
  - `doBreakdown()`: freezes/unfreezes index slices to isolate `PDF_unc` vs (if
    `fitScales`) `scale_unc`, then prints both as separate named terms.
  - `plotFitResults()` has a separate `PDF_postfit` plot (errorbar of post-fit
    value +/- unc per index, `pulls = value/sqrt(1-unc**2)`) -- the analogue the
    user wants for TNP.
- Externalized scale: `doScaleVariations(corr)` -> `getScaleVariationsMeas(meas)` ->
  `setGlobalRescalingFactorTheory(meas, mur, muf)`, which sets
  `self.globalRescalingFactorTheory = readth(mur,muf)/readth(central)` (via
  `self._readth`, self-area-normalized, per-measurement, drops bin 0). `corr=True`
  uses `meas=self.exp_name`, applied uniformly to *all* bins at once -- the "fully
  correlated across measurements" externalized path.
- `th_xsec_nnlo.py`'s `th_xsec.readNormXsecMeas`/`readNormXsec`/`_read_dist`: read
  one measurement's self-area-normalized rho, drop bin 0, concatenate across
  measurements. `_find_energy_obj`/`_find_order_key`/`_find_mass_key`/
  `_find_pdf_family`/`_meas_obs_label`/`_scale_key` are the JSON-traversal helpers
  to reuse.
- Numerically verified (this session, `mtpole-ttj-pyconvino/plotTNPandExpandedScale.py`):
  self-area normalization and xsec-ratio normalization agree to ~8e-6 everywhere
  (floating-point noise) -- `xseccms`/`xsecatl` is exactly the bin-width-integrated
  sum of the rho histogram's own bins in this theory calculation. TNP width=1 is the
  raw stored eigenvector step = 68% CL (matches a unit-Gaussian-prior nuisance
  directly, no extra rescale needed).

## Design

### 1. `th_xsec_nnlo.py`: new reader methods (mirror `_read_dist`/`readNormXsecMeas`/`readNormXsec`)

- `readTNPVariationMeas(PDF, mass, meas, ev_idx, order)` /
  `readTNPVariation(PDF, mass, ev_idx, order)`: read
  `tnp_eigenvectors[family][obs_label][ev_idx]` at the central scale
  (`_scale_key('1.0','1.0')`), plus the central member's absolute
  `distributions[obs_label]` for area-normalization. Symmetric two-sided formula
  (same one used for asymmetric-Hessian PDF pairs, `(S_plus - S_minus)/2`):
  ```
  norm(x) = x / sum(x * binwidth)
  shift = (norm(central + delta_i) - norm(central - delta_i)) / (2 * norm(central))
  ```
  drop bin 0 (matches `_read_dist`), concatenate across `self.exp_input.input_meas`
  for the "Variation" (no-`Meas`) version.
- `numTNPVariations(PDF, mass, order)`: `len(tnp_eigenvectors[family][obs_label])`
  (3/6/9 for LO/NLO/NNLO).
- `readExpandedNormXsecMeas(PDF, mass, meas, mur, muf, order)` /
  `readExpandedNormXsec(...)`: same traversal as `_read_dist` but read
  `normalized_expanded[family][obs_label]['binContent']` directly (already a
  normalized ratio), drop bin 0, concatenate.
- All raise a clear `KeyError` if family/order isn't present (mirrors
  `_find_pdf_family`'s "surface a clear error, don't silently fall back"
  philosophy) -- this is what enforces "PDF4LHC21_40 only, LO/NLO/NNLO only"
  without a separate special-cased check.

### 2. `fit_object.py`: TNP nuisances (`--tnp`), as a fully separate block from PDF

- New constructor params `tnp=False, expScale=False` -> `self.tnp`, `self.expScale`.
- New `_buildTNPVariationShifts()`: calls the new readers, returns `(names,
  shifts)`, `shifts` shape `(n_tnp, n_bins)` -- same fractional-shift convention
  as `self.PDF_vars`, but stored as **`self._tnp_var_names` / `self.TNP_vars`**,
  never touching `self._pdf_var_names`/`self.PDF_vars`.
- Layout becomes `[mass][PDF][TNP if --tnp][scale if --fitScales]` (mutually
  exclusive, so at most one of the last two blocks exists -- see validation
  below). Every place that currently branches on `self.fitScales` to
  append/read/prior/free its own trailing block gets a parallel `if self.tnp:`
  branch reading `self.TNP_vars`/`self._tnp_var_names` instead, at index range
  `[nMassParams+len(PDF_vars), nMassParams+len(PDF_vars)+len(TNP_vars))`:
  - `initialiseMinuit()`: append `zeros(len(TNP_vars))` when `self.tnp`.
  - `_chi2()`: add `sum(params[tnp_lo:tnp_hi]**2)` to `priors` when `self.tnp`
    (unit Gaussian prior, un-rescaled -- TNPs have no CL-rescale concept, unlike
    some PDF families).
  - `getInterpolatedXsec()`'s non-`_mdv` branch: add a
    `TNP_weights = array([self.TNP_vars[i]*params[tnp_lo+i] for i in range(len(self.TNP_vars))])`
    block, multiplied in exactly like `PDF_weights`/`scale_weights` (respecting
    `self.multiplicative`).
  - `doBreakdown()`: **needs one correctness fix** -- its current unconditional
    `minuit.fixed[self.nMassParams+len(self.PDF_vars):] = False` (used twice,
    computing the PDF-frozen baseline for `PDF_unc`) implicitly relies on nothing
    occupying that trailing slice unless `fitScales`. Once TNP can occupy it, this
    must become `if self.fitScales: minuit.fixed[...] = False` explicitly, so TNP
    correctly stays frozen while `PDF_unc` is computed (otherwise `PDF_unc` would
    silently absorb TNP's contribution -- exactly the conflation to avoid). Then
    add a **new**, `scale_unc`-mirroring block: `if self.tnp:` freeze only the TNP
    slice, migrad, `self.tnp_unc = sqrt(self.minuit.errors[mass]**2 -
    tnp_frozen_errors[mass]**2)`, and extend the printed summary line to include
    `+/- {:.2f} (TNP)` when `self.tnp` (parallel to the existing `(scale)` term).
  - `fitMassMinuit()`: skip `doScaleVariations` entirely when `self.tnp` (replaced
    by the TNP nuisance treatment).
  - `_write_results_json()`: write `mass_unc_tnp` analogous to `mass_unc_scale`,
    guarded by `hasattr(self, 'tnp_unc')`.
- `plotFitResults()`: **existing `PDF_postfit` block is untouched** (it already
  only iterates `len(self.PDF_vars)`, which no longer includes TNP). Add a new,
  separate `TNP_postfit` block right after it -- same errorbar+pulls style
  (`pulls = value/sqrt(1-unc**2)`), x-axis "TNP eigenvector (order=...)", saved as
  `TNP_postfit.pdf`, gated on `self.tnp`, indexing `params_w_cov[tnp_lo:tnp_hi]`.
- Not touched (explicitly out of scope, flagged with a code comment where
  relevant): `doPDFEigenvectorBreakdown`, `doSystEigenBreakdown`, `doNuisanceFit`
  stay PDF-only / don't need to know about TNP unless `--tnp` is combined with
  those opt-in flags, which the validation below forbids for now -- avoids
  silently wrong labeling in code paths not exercised by this change.

### 3. `fit_object.py`: expanded scale variation (`--expScale`)

- Add `expanded=False` param to `setGlobalRescalingFactorTheory`; when true,
  source values from `self.th.readExpandedNormXsec(Meas)` instead of
  `self._readth`, at both call sites inside that one function. Thread
  `expanded=self.expScale`. `getScaleVariationsMeas`/`doScaleVariations` need
  **no changes** -- they call `setGlobalRescalingFactorTheory` and stay agnostic
  to which source it used.

### 4. `doFit.py`: CLI + validation

- New flags `--tnp`, `--expScale` (`action='store_true'`).
- Validation (mirrors existing style, e.g. `--fitScales`+`--freezePDFs`):
  - require `--stripper`
  - mutually exclusive with each other (both explicitly "replace the externalized
    scale var")
  - incompatible with `--fitScales` (shares the same trailing param slot
    conceptually), `--simpleFit` (Minuit-only concept), `--massDepVariations` (no
    per-mass TNP/expanded interpolation built -- future extension), and (for now)
    `--pdfEigen`/`--systEigen`/`--nuisanceFit`/`--budget`/`--globalImpacts` for
    `--tnp` specifically (those breakdown paths aren't TNP-aware; forbidding
    avoids silently mislabeled output)
  - `--order` must be LO/NLO/NNLO (not `NNLO_noVVF`, which has no TNP/expanded data)
  - `--PDF` must resolve (via `configs.PDF_dict`) to `PDF4LHC21_40`
- Pass `tnp=args.tnp, expScale=args.expScale` into `fit(...)`.

## Verification

- `python3 -m py_compile` all three changed files.
- Run `doFit.py --stripper --order NNLO --PDF PDF4LHC21_40 --tnp` against the real
  `--exp` NPZ input already used in this repo; confirm: Minuit converges,
  `TNP_postfit.pdf` + unchanged `PDF_postfit.pdf` both written, printed breakdown
  line shows a separate `(TNP)` term (not folded into `(PDF)`), no
  `doScaleVariations` log lines, nuisance count matches NNLO's 9 TNP eigenvectors
  (`--order LO`/`NLO` -> 3/6).
- Compare `self.PDF_unc` from a `--tnp` run against a plain (no `--tnp`, no
  `--fitScales`) run at the same order/PDF -- should be numerically identical
  (proves the `doBreakdown` fix actually stopped TNP from leaking into the
  PDF-frozen baseline).
- Run `doFit.py --stripper --order NNLO --PDF PDF4LHC21_40 --expScale`; confirm
  externalized scale-variation output is produced as before, and diff the printed
  "Correlated scale uncertainty" number against the same run *without*
  `--expScale` -- expect a small (few-%) but non-zero difference (sanity check
  it's reading different data, consistent with the earlier
  normalized_expanded-vs-classic-ratio numeric check).
- Confirm each guarded invalid combination (`--tnp --fitScales`, `--tnp
  --expScale`, `--tnp --order NNLO_noVVF`, `--tnp --PDF CT18NNLO`, `--tnp
  --massDepVariations`, `--tnp --pdfEigen`, any without `--stripper`) raises the
  intended `ValueError` before any fit work starts.

## Progress log

- 2026-08-04: plan approved by user (with one correction from an earlier draft:
  keep TNP as a fully separate parameter block from PDF, not merged into
  `self.PDF_vars`, for clarity). Nothing implemented yet.
- 2026-08-04: Implemented and verified against real `--exp` NPZ input + the
  `ttbarj-nnlo-cms-atlas-analysis/newdata.json` theory JSON (the one with
  `tnp_eigenvectors`/`normalized_expanded` -- NOT `theory-data/newdata.json`,
  which is stale and lacks both fields; pass `--stripperPath` explicitly to
  point at the former until it's copied over). `numTNPVariations` confirmed
  3/6/9 for LO/NLO/NNLO. All invalid-combination guards raise before any fit
  work starts.
  **Correction to this plan's `doBreakdown()` fix, found by actually running the
  numbers (Verification step 3, "PDF_unc should be numerically identical
  between `--tnp` and a plain run"):** the design above ("must become `if
  self.fitScales: minuit.fixed[...] = False` explicitly, so TNP stays frozen")
  was implemented as written and is WRONG -- it produces `PDF_unc` that's
  nearly the quadrature sum of the plain run's `PDF_unc` and `tnp_unc`
  (0.475 vs plain 0.244, tnp_unc 0.409, sqrt(0.244^2+0.409^2)=0.476), i.e. it
  causes exactly the PDF/TNP conflation it claimed to prevent. The correct fix
  is to leave the trailing-slice unfreeze **unconditional** (i.e. don't touch
  this code at all -- the original pre-plan code was already correct): PDF_unc
  isolates PDF's own marginal contribution by comparing "PDF frozen" vs "PDF
  free" while every other block (scale if `--fitScales`, TNP if `--tnp`) keeps
  the *same* free/frozen status in both the baseline and the full fit --
  exactly mirroring how `scale_unc`/`tnp_unc` isolate their own block by
  freezing only themselves. With the unconditional version, `--tnp`'s PDF_unc
  (0.208) is close to the plain run's (0.244) -- not exactly equal (real
  mass-PDF-TNP correlation), but consistent with the pre-existing, accepted
  imprecision of this quadrature-subtraction breakdown (verified: `--fitScales`
  alone already doesn't close to `mass_unc_total` in quadrature either, ~6%
  off, same as this ~15%). Implemented in `fit_object.py`'s `doBreakdown()` as
  the original unconditional `minuit.fixed[self.nMassParams+len(self.PDF_vars):]
  = False` (both occurrences), now with a comment explaining why it must stay
  unconditional -- do not reintroduce the `if self.fitScales:` guard.

- 2026-08-04 (review pass, same day): four corrections/extensions, all driven by
  real fit numbers rather than reasoning.

  **(a) The `doBreakdown()` conclusion above was half right and half wrong.** Keeping
  the trailing-slice unfreeze unconditional is correct *for `PDF_unc`* (it is what makes
  PDF_unc a marginal term, comparable to a plain run). It is wrong for `exp_unc`/
  `interp_unc`, which used the SAME baseline and therefore silently absorbed the whole
  TNP contribution -- which `tnp_unc` then reported a second time. Measured on
  `--exp ATLAS_813TeV_CMS_13TeV_npz --order NNLO --PDF PDF4LHC21_40 --polyOrder 1 --tnp`:
  exp 0.974 (+) interp 0.137 (+) PDF 0.241 (+) TNP 0.465 = 1.114 in quadrature against a
  HESSE total of 1.013, a 10% overcount; `doExpStatSystBreakdown`'s stat/syst split, being
  anchored to `exp_unc`, inherited it (0.488/0.843 instead of 0.427/0.736). The same
  defect pre-existed for `--fitScales`. Fixed: `PDF_unc`/`tnp_unc`/`scale_unc` keep their
  marginal (freeze-only-my-own-block) definitions, while the exp/interp baseline now
  freezes EVERY nuisance. After the fix: exp 0.851 / interp 0.116 / PDF 0.241 / TNP 0.465
  = 1.006 vs 1.013 (-0.68%). The plain no-TNP/no-fitScales path is bit-identical to
  before (the trailing slice is empty there -- which is exactly why the old code looked
  right). The remaining -0.68% is a genuine PDF<->TNP overlap and is now printed as a
  `breakdown closure:` line + `mass_unc_breakdown_{quadrature_sum,closure_residual}` in
  results.json, instead of being invisible.

  **(b) Marginal terms cannot close in general**, so a new `doAdditiveMassBreakdown()`
  gives the decomposition that does. Linear/Gaussian global impacts (arXiv:2307.04007) on
  the theory side, mirroring `doGlobalImpacts` on the experimental side:
  `sigma_m^2 = w.Cexp.w + w.Cinterp.w + sum_k (w.d_k)^2`, `w = Ceff^-1 J/(J.Ceff^-1.J)`,
  `Ceff = Cexp + Cinterp + sum_k d_k d_k^T`. `J` and `d_k` are central differences of
  `getInterpolatedXsec` at the post-fit point, so the multiplicative composition, the PDF
  priors' CL rescale and `--massDepVariations` are all included exactly as the fit uses
  them; a single `sigma_HESSE/sigma_linear` anchor makes it close to the fit's own total
  while preserving additivity. Verified to close to 6 digits: plain 0.8565/0.1066/0.2351
  -> 0.894578 = HESSE; `--tnp` 0.8517/0.0873/0.2351/0.4867 -> 1.012524 = HESSE. Outputs:
  `mass_breakdown_additive.csv` (block totals + per-nuisance impacts), a per-nuisance
  bar plot, and `mass_breakdown_additive` in results.json. Runs unconditionally, cheap.

  **Closure is STRUCTURAL** (`total == anchor * sigma_lin == sigma_HESSE` identically), so
  it holds for every PDF set, poly order and flag combination -- confirmed on CT18NNLO
  (1.02736580 vs HESSE 1.0273658317, 8 digits). What varies is the ANCHOR, i.e. how far
  the fit sits from the linear limit. Measured 2026-08-05 by numerically differentiating
  the full chi2 Hessian at the post-fit point (PDF4LHC21_40, NNLO, ATLAS_813TeV_CMS_13TeV_npz):

  | poly | tnp | mdv | sigma_HESSE | sigma_lin | anchor | GN vs lin | max abs dH mass rows |
  |------|-----|-----|-------------|-----------|--------|-----------|----------------------|
  | 1 | -   | -   | 0.894578 | 0.904120 | 0.9894 | 1.000000 | 0.0328 |
  | 1 | yes | -   | 1.012523 | 1.052667 | 0.9619 | 1.000000 | 0.0441 |
  | 1 | yes | yes | 1.056828 | 1.053230 | 1.0034 | 1.000000 | 0.0186 |
  | 2 | -   | -   | 0.850049 | 0.837160 | 1.0154 | 1.000000 | 0.0692 |
  | 2 | yes | -   | 0.941104 | 0.951428 | 0.9891 | 1.000000 | 0.0556 |
  | 2 | yes | yes | 0.977647 | 0.953937 | 1.0249 | 1.000000 | 0.0550 |

  `GN vs lin = 1.000000` everywhere: the `Ceff` formula reproduces the Gauss-Newton Hessian
  EXACTLY, so the anchor is never algebra error -- it is entirely the d2model terms that
  Gauss-Newton drops (`sigma_exact` from the numeric Hessian matches HESSE to <=1e-5 in all
  six). Sources: the multiplicative composition `xsec *= prod(1 + w_i theta_i)` (separately
  for PDF and TNP, so the two blocks also multiply each other) and the mass polynomial.
  Going poly1 -> poly2 roughly DOUBLES the mass-row second derivatives (the new d2x/dm2)
  and flips the anchor's sign. `--massDepVariations` shrinks the PDF-row terms (0.0299 ->
  0.0137 at poly2) but not the mass-row ones, so it helps a lot at poly1 (3.8% -> 0.34%)
  and HURTS at poly2 (1.09% -> 2.49%) -- the flags interact; there is no single best
  setting. Note also sigma_HESSE itself moves ~7% between poly1 and poly2 (1.0125 vs
  0.9411 with `--tnp`), larger than any breakdown effect here; `--polyOrder` defaults to 2
  while this whole study used 1.

  **(c) `meta-data/setups.py` was never needed here, and the restriction was on THREE
  axes, not one.** `convertData.py::add_cross_order_quantities` drove its loops from
  setups.py, which enumerates only the `PDF4LHC21` group's three (HT4/HT2/HT) scale
  choices -- but setups.py supplied only two things: the `scale_key` string (available
  directly from `mass_entry['scales']`) and the `CC/UU/DD/...` label (used solely to gate
  TNPs to the CC points). Both are redundant, so setups.py has been dropped from this
  script entirely, along with the `setup_info_by_mass` plumbing that fed it. The
  dependency was inherited from `compute-tnp-eigenvectors.py`, which hardcodes `_CC` in
  its `data/*.dat` filenames and loops `for pdf in setups.keys()` -- a limitation of what
  was extracted into those .dat files, NOT of the TNP model, which needs nothing beyond
  one family's own central predictions for the contributions in its order chain, at one
  scale point. What this unlocked, all verified present in the XMLs:
    - **families**: 1 -> 11 (every family has central members at these scale points);
    - **scale points**: TNPs were written at 3 of 13, now all 13
      (`normalized_expanded` already covered all 13 -- its loop had no CC gate);
    - **orders**: `NNLO_noVVF` was skipped only because the hardcoded
      `CONTRIBUTIONS = ['LO','NLO','NNLO']` list drove the loop and gave it no chain.
      Its series is well defined as `[LO, NLO-LO, NNLO_noVVF-NLO]`, and since
      NNLO_noVVF is the *complete physical prediction* while NNLO is the bare
      epsilon-pole expansion, that is arguably the more meaningful TNP of the two.
  Index arithmetic (`CONTRIBUTIONS.index(order)`) replaced by an explicit `ORDER_CHAINS`
  table; `compute_tnp_eigenvectors`/`compute_normalized_expanded` take the chain.
  Result: 10440 -> 19920 derived entries (exactly 2490 per order per quantity =
  83 (scale, family) combos x 3 obs x 5 masses x 2 energies), 189.1 -> 217.8 MB, with
  **zero pre-existing entries lost and max abs diff 0.0** on all of them.
  **The regenerated JSON must use `--smeared`** -- that is the production convention and
  the default is unsmeared; regenerating without it shifts every distribution by 0.3-2%
  and every TNP by up to 9%. `theory-data/newdata.json` has been replaced with the new
  file (verified a bit-identical superset of the old one: 115440/115440 common entries
  equal), so the `--stripperPath` gotcha above no longer applies. `doFit.py` now only
  rejects legacy-ROOT-only `--PDF` values; the `--order in (LO,NLO,NNLO)` guard is gone.

  **(d) `--massDepVariations` now works with `--tnp` and `--expScale`.** The theory JSON
  carries `tnp_eigenvectors` at all 5 mass points, so there was never a data reason to
  freeze the TNP fractional shift at `default_mass`. `th_xsec_nnlo` gained
  `readTNPVariationValueWithError` (pseudo-absolute `central + shift`, same convention as
  `_readPDFVariationValueWithError`), `_interpolateVariations` fits it per bin across mass
  into `self.xsec_interp_tnpvar`, and `getInterpolatedXsec`'s `_mdv` branch uses it. Under
  `--expScale` the scale polynomials are now built from the expanded ratio too, otherwise
  the mass-dependent rescale would mix an expanded central ratio with classic per-mass
  scale shapes. Both incompatibility guards dropped from `doFit.py`.

  **(e) Chebychev0 model parameters 3, 6, ... are exact null directions.** `tnps.py`'s
  `get_prediction_next_term` subtracts `theta[it*nppt]` straight back off the Chebyshev
  sum for every order it>=1, so each correction block's constant coefficient cancels
  analytically. NNLO stores 9 eigenvectors but only 7 can move the prediction (NLO 6->5,
  LO 3->3); the other two enter Minuit as parameters pinned by their prior alone.

  This is a property of the MODEL, not of these inputs, and the null parameters are
  **kept** (user decision, 2026-08-05): `_chi2` gives every TNP a unit-Gaussian prior, so a
  null direction simply sits at 0 with error 1 and contributes exactly 0 to every
  breakdown, while the nuisance index stays equal to the model-parameter index with no
  gaps at any order. Verified a numerical no-op: keeping all 9 reproduces the 7-parameter
  run to the last printed digit (mt 171.709096 +/- 1.012524; exp 0.8517 / interp 0.0873 /
  PDF 0.2351 / TNP 0.4867; anchor 0.9619), and `mass_breakdown_additive.csv` lists
  `tnp_NNLO_3` and `tnp_NNLO_6` with impact exactly 0. Two consequences handled:
  `plotPostFitParams`'s pull denominator `sqrt(1 - unc^2)` vanishes for an exactly
  unconstrained nuisance, so the pull is reported as 0 there; and
  `plotTNPandExpandedScale.py` draws the null variations (labelled `(null)`, flat on
  central) instead of hiding them.

  Note this also means `tnps.py`'s `model_n_parameters` = `(sigma_maxpower+1)*(poly+1)`
  over-counts the genuinely free directions by `sigma_maxpower`; the true count is
  `sigma_maxpower*poly + poly + 1` (9->7, 6->5, 3->3).

  **Open question, not a bug:** the TNP eigenvector basis is defined in BIN-INDEX space
  (`model_x` = midpoints of `linspace(-1, 1, nbins+1)`), and rhocms13/rhoatlas13/rhoatlas8
  have 4/7/8 bins. Tying eigenvector i across all three measurements and both energies --
  which is the intended "fully correlated" design -- therefore correlates the same
  Chebyshev coefficient in each observable's own bin-index variable, not literally the
  same function of rho. Defensible, but worth deciding deliberately whether the
  correlation should instead be per-observable or per-energy.
