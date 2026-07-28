# Fit-matrix status — 2026-07-13 (updated, fourth session same day)

Written so a clean session can resume without re-deriving today's work. Full
detail lives in project memory (`fit_matrix_storage_blinding_plan`,
`condor_fitmatrix_eigh_crash`, `condor_fitmatrix_splitmasses_crash`,
`condor_theory_data_afs_load`, `afs_environment_gotchas`,
`feedback_verify_all_axes_plotted`) — this file is a shorter pointer.

## Fourth-session update (same day)

- Phase 4c row-mapping confirmed by user; `loadFromMatrix.py` now maps all 4
  combination rows (was 3). New `MtopSummaryPlot/matplotlib_paperplot.py`
  (matplotlib/mplhep, no ROOT/cppyy) renders the paper summary plot with
  these 4 real rows, `makePlot.py`/`measurementData.py` left untouched.
- New uncertainty-breakdown plot in `matrix_comparison_plot.py`
  (`plot_unc_breakdown`): exp/+PDF/+interp single markers per point
  (quadrature-cumulative). 4th "+scale" marker has no data yet — see below.
- `fit_object.py` fix: the externalized correlated scale uncertainty
  (`doScaleVariations(corr=True)`, always computed, NOT `--fitScales`) is
  now persisted to `results.json` (`mass_unc_scale_ext_up/down`) instead of
  print-only. **Only affects future fits** — the existing 232-job matrix
  predates this fix and needs a partial rerun (224 single-mass jobs) to
  populate it. Not yet triggered — needs user go-ahead.
- New `mtpole-ttj-pyconvino/matrix_2d_plot.py`: matrix-based 2-POI/3-POI
  ellipse plots + single-vs-split Δχ² compatibility p-value, for all 4
  POI-split entries in the matrix (supersedes old `2D_plot.py` for this,
  which silently dropped the 3-POI case and had a dataset_key grouping bug).
- All of the above verified against the real `fitmatrix_results_20260713/`
  matrix and visually reviewed. Still nothing committed.

## Third-session update (same day)

- Real bug the user caught by inspection: every diagnostic plot silently
  only showed `poly_order=2`, even though the whole 232-job matrix
  deliberately ran BOTH poly1 and poly2 for every axis-tuple. Fixed in
  `matrix_comparison_plot.py`: `main()` now loops both poly orders, writing
  each into its own `<dataset_key>/poly{1,2}/` subfolder, plus a new
  `<dataset_key>/poly_comparison.pdf` per fit (all 16 order×PDF points,
  polyOrder=1 vs 2 shown side by side on the same rows). Regenerated: 49
  plot files total (7 fits × (2 poly × 3 plots + 1 poly_comparison)),
  spot-checked visually. See `feedback_verify_all_axes_plotted` memory for
  the general lesson.
- Cleaned up `condorJobs_20260713_fitmatrix_v3/` (223MB local scaffolding
  for the now-completed/collected 232-job run — results safe in
  `fitmatrix_results_20260713/matrix.csv` + EOS, 2.7GB confirmed still
  there). AFS quota back to 78%.

## Second-session update (same day)

- Matrix extended from 1 fit (full combo only) to all **7 user-facing
  fits** (ATLAS8-only, ATLAS13-only, CMS-only, ATLAS8+13, ATLAS8+CMS13,
  ATLAS13+CMS13, full combo) each getting the full stripper order×PDF
  sweep. **Real 232-job Condor run: 232/232 passed**, job dir
  `condorJobs_20260713_fitmatrix_v3/`, EOS base `.../condor/_20260713_fitmatrix_v3/`.
- Item 2 below (missing `ATLAS_813TeV_npz`/single legacy point) is now
  fixed — new `atlas_energy_combo` category in `generateFitMatrix.py`.
- New per-fit diagnostic plots: `mtpole-ttj-pyconvino/matrix_comparison_plot.py`
  now loops all 7 fits, plus a new `legacy_vs_new_nlo` plot type (both a
  dedicated plot and an extra row in the PDF-to-PDF NLO panel). Results:
  `fitmatrix_results_20260713/comparison_plots/<dataset_key>/{order_to_order,pdf_to_pdf,legacy_vs_new_nlo}.{pdf,png}`.
- New `MtopSummaryPlot/makePaperPlot.py` (Version A/B paper plot, separate
  from `makePlot.py` which stays untouched) — row-selection logic verified
  against real data, but the actual ROOT-rendered image is **not yet seen**
  (this shell's PyROOT is missing `cppyy`, a pre-existing env gap — run it
  yourself wherever `makePlot.py` normally works).
- Found+fixed a real tarball-duplication bug in `generateFitMatrix.py`
  (was rebuilding the same ~31MB tarball ~30x per dataset) — now dedups via
  symlink + `--reuse-tarball` (232-job local footprint: 222MB, not ~7GB).
  AFS does not support cross-directory hardlinks — see `afs_environment_gotchas`.
- Cleaned up: old `condorJobs_20260712_fitmatrix_v2/` local scaffolding
  (1.3GB, results already safe elsewhere) and the failed/partial v3 attempts.

## What's done (original, first session)

- **Real 50-job Condor run succeeded: 50/50 passed.** Job dir
  `condorJobs_20260712_fitmatrix_v2/`, EOS base
  `/eos/cms/store/group/cmst3/group/top/sewuchte/MassCombination/condor/_20260712_fitmatrix_v2/`.
  25 axis-tuples × `--polyOrder 1`/`2` (new axis, added this session).
- Two real bugs found and fixed on this run (both root-caused and verified,
  not just theorized — see the two `condor_fitmatrix_*_crash` memories):
  1. `eigh`/`LinAlgError` crash (was blocking the *previous* session's
     attempt) — `--do-impacts` now defaults on in `generateFitMatrix.py`.
  2. `--splitMasses` `IndexError` (new this session, 6/50 jobs) — one-line
     `np.asarray()` fix in `fit_object.py`'s
     `get_correlation_matrix_from_covariance_matrix`.
- NNLO/`--stripper` theory data (`theory-data/newdata.json`, 115MB) now also
  ships in the per-job tarball instead of reading live off AFS — mirrors the
  NLO leg that was already done. `PLAN_afs_load_fix.md` (the old design doc)
  was deleted since both legs are now complete.
- `tools/collectFitMatrix.py` and `mtpole-ttj-pyconvino/matrix_comparison_plot.py`
  ran against real data (not synthetic fixtures) for the first time — found
  + fixed a `pd.read_csv` dtype bug along the way (legacy numeric PDF IDs
  silently become `int64` on CSV round-trip).
- New `MtopSummaryPlot/loadFromMatrix.py` (Phase 4c): turns the aggregated
  matrix into real, salt-blinded rows for 3 of the 4 "boosted"/combination
  rows in `measurementData.py`.
- Results: `fitmatrix_results_20260712/` — `matrix.csv`, `matrix_meta.json`,
  `comparison_plots/{order_to_order,pdf_to_pdf}.{pdf,png}`,
  `measurementData_matrix.py`.
- Cleanup: original failed run's evidence deleted (local + EOS, 841M + 323M).

## What's NOT done — pick up here

1. **Phase 4c is not wired into the live plot.** `loadFromMatrix.py`'s
   output (`fitmatrix_results_20260712/measurementData_matrix.py`) maps 3
   matrix axis-tuples to `measurementData.py`'s hand-typed "boosted" rows
   (`"LHC "`, `"ATLAS + CMS"`/`"13"`, `"ATLAS + CMS"`/`"8+13"`) — this
   mapping is the agent's inference from terse labels and was **deliberately
   not applied** to the live, paper-facing `measurementData.py`/`makePlot.py`
   without user confirmation first. Confirm the mapping, then decide how to
   wire it in (`loadFromMatrix.py` writes directly into
   `measurementData.py`, or `measurementData.py` imports
   `matrix_measurements` from the generated file — not yet decided).
2. ~~4th "boosted" row (`"ATLAS 8+13"`) has no matching matrix entry~~ —
   **FIXED second session**: `atlas_energy_combo` category added and the
   232-job run includes it. `loadFromMatrix.py` itself hasn't been re-run
   against the new matrix yet to pick this 4th row up, though (still only
   maps 3 of 4 as written).
3. **Nothing from any of the three 2026-07-13 sessions is committed.**
   Touches (on top of the first-session list): `mtpole-ttj-pyconvino/matrix_axes.py`
   (`SWEEP_DATASET_KEYS`, `DATASET_DISPLAY`), `condorManagement/generateFitMatrix.py`
   (7-way sweep, tarball dedup), `mtpole-ttj-pyconvino/matrix_comparison_plot.py`
   (per-fit loop, `plot_legacy_vs_new`, poly1/poly2 loop + `plot_poly_comparison`),
   `MtopSummaryPlot/makePaperPlot.py` (new). Decide commit strategy (one
   combined commit vs. topically split).
4. Blind salt (`~/.masscomb_blind_salt`) confirmed intact and correct
   (fingerprint `4de0c068e866b4fb`) — **not** stored inside the repo
   (deliberate, checked explicitly this session).

Once this file's content is stale (Phase 4c wired in + committed), delete it
— it's a handoff note, not a permanent doc.
