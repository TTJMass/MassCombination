# TTJetMassCombination

Extraction of the top quark pole mass (m_t) by combining differential
ttbar cross-section measurements from ATLAS and CMS, then fitting NLO/NNLO
theory predictions. This repo centralizes the combination setup, fitting,
plotting scripts, and Condor batch infrastructure.

## Submodules

| Directory | Purpose |
|---|---|
| [`pyconvino/`](pyconvino/README.md) | Python/JAX statistical combination (χ² built with JAX, minimized with SciPy). Produces `*_result.{txt,npz,json}`. |
| [`mtpole-ttj-pyconvino/`](mtpole-ttj-pyconvino/README.md) | Top-mass extraction: fits pyconvino's combination output against NLO/NNLO theory predictions. |
| [`MtopSummaryPlot/`](MtopSummaryPlot/README.md) | Matplotlib/mplhep summary ("forest") plot of top mass measurements. |
| [`theory-data/`](theory-data/README.md) | NLO/NNLO POWHEG theory predictions (8/13 TeV). |

Note: `mtpole-ttj-pyconvino` shares its remote with `mtpole-ttj` (an older
branch of the same repo); it's checked out as a second submodule on its own
branch, not a separate GitHub project.

## How to install

```bash
git clone --recursive git@github.com:TTJMass/MassCombination.git
```
or
```bash
git clone git@github.com:TTJMass/MassCombination.git
cd MassCombination
git submodule update --init --recursive
```

Each submodule has its own Python dependencies — see its README for the
`pip install` command. All of them run inside the shared `masscomb` conda
env for Condor batch jobs (see `condorManagement/rebuild_env_pack.sh`).

## Inputs

- **Experimental**: `Inputs_ATLAS/` and `Inputs_CMS/` hold the reference
  experimental measurement/correlation files; the per-setup combination
  configs under `pyconvino/ConvinoSetups/*/` are built from these.
- **Theory**: NLO/NNLO POWHEG predictions live in the `theory-data`
  submodule — see [`theory-data/README.md`](theory-data/README.md) for how
  to update them.

## How to run

1. **Combine**: run `pyconvino` on a `ConvinoSetups/<setup>/rho_config.txt`
   to produce a combined result (see `pyconvino/README.md`).
2. **Fit**: run `mtpole-ttj-pyconvino/doFit.py` on the combination output
   against theory to extract m_t (see `mtpole-ttj-pyconvino/README.md`).
3. **Plot**: run `MtopSummaryPlot/matplotlib_paperplot.py` to render the
   summary forest plot.

For batch running on HTCondor, `condorManagement/createCondorJobs.py`
packages a setup + the pyconvino/mtpole-ttj-pyconvino stack (`--mode new`,
the default) and submits scan jobs; `writeAllCondorJobs_pyconvino.sh` drives
it across all `ConvinoSetups`. `condorManagement/HTC_monitor.py` monitors
and resubmits failed jobs. `tools/collectCondorScans.py` and
`tools/collectFitMatrix.py` aggregate the results afterwards.

## docs/

Background/derivation notes (e.g. `docs/GLOBAL_IMPACTS_MATH.md` on the
additive stat/syst/per-group mass-uncertainty breakdown) live in `docs/`.

## Decommissioned old stack

The original C++/ROOT combination (`Convino/`), the ROOT-based fit
(`mtpole-ttj/`), and `MtopSummaryPlot`'s PyROOT plotting path were removed
from git on 2026-07-28 in favor of the pyconvino/mtpole-ttj-pyconvino/
matplotlib stack above. Local tarball archives of the removed content are
kept in `archive/` (gitignored, not part of the repo).
