# Unblinding page (8 Oct 2026)

Click-to-reveal page of the unblinded results (the main result first).

- `build_page.py blind|real frozen|lep075 OUT`: reads the fit matrix, removes the blinding salt
  (`~/Claude_MassComb/.masscomb_blind_salt`, fingerprint 4de0c068e866b4fb) and writes the values and
  plots straight into `OUT/unblinding.html`. Prints no value. `blind` keeps the salt (for checks).
  `lep075` reads `matrix_lep075.csv` and `scanrun_lep075/` from `DATA` (set in the script).
- `capture_scan.py`: runs `doFit.py` unchanged (blinded) and stores the chi2 scan points
  (`scan_capture.json`), used for the chi2 scan plot.
- `template.html`: page layout.
- The summary plot comes from `MtopSummaryPlot/prl_summary_options4.py`.

Plot scripts remove the salt at load time with `blinding.unblind_matrix` (mtpole-ttj-pyconvino).
Set `MASSCOMB_UNBLIND_DRYRUN=1` to keep the salt (checks without real values).
