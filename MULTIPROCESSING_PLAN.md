# Speeding up `mtpole-ttj-pyconvino` fits via multiprocessing — plan

Status: **analysis done, nothing implemented yet**. Written 2026-08-10.
Resume cold from this file — it's the durable copy of a design discussion that
happened in a Claude Code session (not synced anywhere else).

## Context

Fits take noticeably long for the heavy 3-input combination (~600 PDF/syst
nuisance parameters). We already added `--bindropCheck` (see the `doBindropCheck`
diagnostic in `fit_object.py`) which does its own scan of independent refits, and
in discussing how to speed *that* up specifically, we did a broader pass over the
whole `mtpole-ttj-pyconvino/fit_object.py` fit pipeline looking for every place that
does "the same independent piece of work" multiple times in a loop — since those are
the safe, high-value multiprocessing targets. This file records that survey so the
work can be picked up without re-deriving it.

All line numbers below are as of 2026-08-10 (after the `--bindropCheck` commit) and
will drift — treat them as "grep for the function name to relocate," not gospel.

## Core safety fact, applies everywhere in this file

`fit.fullChi2`/`fit._chi2` (the closure `iminuit.Minuit` minimizes) reads *mutable
instance state* (`self.xsec_interp`, `self.cov_PDF`, `self._cov_inv`,
`self._dropBinIdx`, ...), and several of the loops below temporarily mutate that
state per-iteration (e.g. `doBreakdown`'s no-interp branch zeroes
`self.cov_interp` for the duration of one refit, then restores it).

- **Threading is unsafe** for any of this: concurrent threads share the same
  `self`, so one iteration's temporary mutation would corrupt another's chi2
  evaluation mid-flight. This isn't just a perf argument, it'd silently produce
  wrong numbers.
- **Fork-based multiprocessing is naturally safe**: `multiprocessing.get_context('fork')`
  gives each worker a private copy-on-write clone of the whole process, including
  `self`. A worker can mutate its own copy of `self.cov_interp` freely — it can
  never leak into siblings or the parent. This also means the current manual
  "save/restore" dance some of these functions do becomes unnecessary if the work
  moves into a forked child: the parent's `self` was never touched in the first
  place.
- Caveat: CPython's reference counting touches (dirties) the memory page of any
  object just by *reading* it, so COW savings are real but smaller than the naive
  "big shared read-only blob costs nothing" intuition — still much better than
  each worker re-loading/re-parsing the ~200 MB theory JSON from disk, which is
  the single biggest fixed cost per process (see prior session's research).
- Cap worker count (e.g. `os.cpu_count()`) and pin BLAS/OpenMP thread count to 1
  inside each worker (`OMP_NUM_THREADS=1` etc.) to avoid N processes × M internal
  threads oversubscribing cores.

## Tier 1 — runs on every single fit, biggest lever

These all share one shape: *deepcopy(`self.minuit`) → mutate a fixed setting →
`migrad()` → read result*, already structured as independent work by the existing
code (it already deep-copies before each branch), just executed sequentially.

- **`doBreakdown()`** (`fit_object.py:1120`) — unconditional, runs on every
  `doFit.py` invocation. Sequentially does 2–5 full re-minimizations over all
  nuisances: PDF-frozen (`minuit_nopdf`, ~line 1136), exp+interp (`minuit`, ~line
  1146), no-interp if `cov_interp` is nonzero (`minuit_nointerp`, ~line 1157,
  temporarily zeroes `self.cov_interp`), no-scale if `--fitScales`
  (`minuit_noscale`, ~line 1169), no-TNP if `--tnp` (`minuit_notnp`, ~line 1178).
  Each branch only reads `self` as a starting template — independent otherwise.
  **Highest priority**: unconditional cost on every fit, 5 sequential ~600-param
  migrad calls in the worst case.

- **`doChi2Scan()`** (`fit_object.py:1040`) — also unconditional, once per mass
  parameter. Loops `scan_range` (~40 points, ±2σ in 0.1 GeV steps), deep-copies
  `self.minuit` once (line ~1047), fixes the mass, and calls `migrad()` per point
  (line ~1052). The cleanest possible target — no per-branch state mutation at
  all, just "fixed mass → migrad → read chi2".

- **`doBindropCheck()`** (`fit_object.py:1832`, added this session) — ~10–20
  independent `initialiseFit(dropBinIdx=...) + migrad/hesse/minos` refits. Already
  discussed with the user as a multiprocessing candidate before this broader
  survey; naturally becomes the third user of the same shared helper.

**Recommended approach**: write one reusable helper, e.g.
`_run_parallel_refits(self, jobs)` in `fit_object.py`, built on
`multiprocessing.get_context('fork').Pool`, where each `job` is a small closure/
callable describing "what to mutate + which minuit to use + what to run", and use
it in all three places above (refactoring `doBindropCheck`'s current sequential
loop onto it too instead of leaving it as a one-off).

## Tier 2 — `--massDepVariations` (flag-gated, lower priority)

- **`_interpolateVariations()`** (`fit_object.py:633`) feeding
  **`_fitPolyAcrossMass()`** (`fit_object.py:705`) — loops over every PDF
  eigenvector (up to ~300 for a full replica family) plus the 6 non-central scale
  points plus TNP eigenvectors, each doing a few theory-JSON dict reads and small
  per-bin `np.polyfit` calls.
  - Read-only w.r.t. `self` (no mutation inside the loop body) — **safe under
    threads**, unlike Tier 1.
  - But each item is cheap (dict lookups + tiny linear algebra, no migrad), so
    process-pool spawn/pickling overhead may exceed the savings per item.
  - **Before committing to multiprocessing here: time it.** Try
    `concurrent.futures.ThreadPoolExecutor` first (numpy's LAPACK least-squares
    call does release the GIL for the underlying `dgelsd`/`dgelsy` call), and only
    fall back to fork-based multiprocessing if that doesn't help.
  - Only exercised under `--massDepVariations`, so lower priority than Tier 1
    despite being the flag the user originally asked about.

## Not good multiprocessing candidates

- **`doNuisanceFit()`** (`fit_object.py:3574`) and **`doSystEigenBreakdown()`**
  (`fit_object.py:3190`): each is *one* large nonlinear `migrad()` over
  ~600–900 params (not a loop of independent smaller fits) — this is where the
  "may take several minutes" MINOS cost lives. Can't be split without changing
  the optimizer itself; out of scope for this plan.
- **`doPDFEigenvectorBreakdown()`** (`fit_object.py:3119`): pure linear algebra
  on an already-computed Hessian covariance matrix — already fast.
- **`fitMass()`** (`fit_object.py:2866`, the `--simpleFit` legacy path): loops
  800 mass-grid points, each doing one small `np.linalg.solve`
  (`computeChi2`, `fit_object.py:514`). Cheap enough per-call that this is a
  **vectorization** candidate (batch the 800 linear solves via numpy broadcasting)
  rather than a parallelization one. Legacy/rarely-used path — low priority
  either way.
- **`doUncertaintyBudget()`** (`fit_object.py:1940`) / `doGlobalImpacts()`
  (`fit_object.py:2294`): matrix-vector products over pre-loaded NPZ arrays via
  `_fisher_weights()` — no refit loop, already fast.

## Already handled, no action needed

- **`run_all.py`** already has `--jobs N` (`run_all.py:35`, dispatch at
  `run_all.py:92-104` via `concurrent.futures.ProcessPoolExecutor`) launching
  independent `doFit.py` subprocesses for sweeps — real OS processes, zero
  shared-state risk, already the right design for that layer.
- Condor scans (`condorManagement/generateFitMatrix.py`) are already externally
  parallel across worker nodes.

## Next steps

1. Decide starting point: Tier 1 shared helper (biggest win, touches the core fit
   path used by every invocation — higher review scrutiny warranted) vs. Tier 2
   (`--massDepVariations`, smaller blast radius since it's flag-gated, but needs
   timing data before choosing threads vs. processes).
2. If Tier 1: implement `_run_parallel_refits`, apply to `doBreakdown` and
   `doChi2Scan` first (both unconditional — easiest to validate: compare
   before/after numbers exactly, they should be bit-identical since it's the same
   computation just reordered), then port `doBindropCheck` onto it.
3. Smoke-test on the real 3-input combination fit (the one motivating this) and
   record the wall-clock before/after in this file or a follow-up.
