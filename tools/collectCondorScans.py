#!/usr/bin/env python3
"""Collect condor scan results (fit.log) and plot scan summaries.

Searches under a jobs folder (and optionally an EOS output base) for
`fit.log` files produced by the Condor job wrappers. Parses lines of the
form:

  raw mt = 171.918555 +/- 0.612156 (exp) +/- 0.434377 (PDF) GeV

and builds a dictionary of correlation scans with values and results.

Produces per-scan plots (3 panels) and two 2D summary plots, and
saves a JSON summary.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
from collections import defaultdict
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
hep.style.use(hep.style.CMS)


RAW_MT_RE = re.compile(
    r"raw\s+mt\s*=\s*([0-9.+-eE]+)\s*\+/-\s*([0-9.+-eE]+)\s*\(exp\)\s*\+/-\s*([0-9.+-eE]+)\s*\(PDF\)",
    flags=re.IGNORECASE,
)


def find_fit_logs(base_paths: List[str]) -> List[str]:
    """Find all `fit.log` files under the given base paths.

    Returns absolute paths.
    """
    files = []
    for bp in base_paths:
        print ('Searching for fit.log under', bp)
        if not bp:
            continue
        if os.path.isfile(bp) and os.path.basename(bp) == 'fit.log':
            files.append(os.path.abspath(bp))
            continue
        for path in glob.glob(os.path.join(bp, '**', 'fit.log'), recursive=True):
            files.append(os.path.abspath(path))
    return sorted(list(set(files)))

def find_convino_logs(base_paths: List[str]) -> List[str]:
    """Find all `convino_run.log` files under the given base paths.

    Returns absolute paths.
    """
    files = []
    for bp in base_paths:
        print ('Searching for convino_run.log under', bp)
        if not bp:
            continue
        if os.path.isfile(bp) and os.path.basename(bp) == 'convino_run.log':
            files.append(os.path.abspath(bp))
            continue
        for path in glob.glob(os.path.join(bp, '**', 'convino_run.log'), recursive=True):
            files.append(os.path.abspath(path))
    return sorted(list(set(files)))


def parse_fit_log(path: str) -> Optional[dict]:
    """Parse a `fit.log` file and extract central and uncertainties.

    Returns dict with keys: central, exp_unc, pdf_unc, total_unc
    or None if not found.
    """
    try:
        txt = open(path, 'r', errors='ignore').read()
    except Exception:
        return None

    m = RAW_MT_RE.search(txt)
    if not m:
        return None

    central = float(m.group(1))
    exp_unc = float(m.group(2))
    pdf_unc = float(m.group(3))
    total_unc = math.hypot(exp_unc, pdf_unc)

    return {
        'central': central,
        'exp_unc': exp_unc,
        'pdf_unc': pdf_unc,
        'total_unc': total_unc,
        'path': os.path.abspath(path),
    }


def _find_sublist_index(hay: List[str], needle: List[str]) -> int:
    """Return index where `needle` sequence appears in `hay`, or -1."""
    if not needle:
        return -1
    for i in range(len(hay) - len(needle) + 1):
        ok = True
        for j in range(len(needle)):
            if hay[i + j] != needle[j]:
                ok = False
                break
        if ok:
            return i
    return -1


def infer_scan_and_value_from_path(path: str, jobs_folder: Optional[str], eos_output: Optional[str] = None) -> (str, Optional[float]):
    """Infer scan name and correlation value from the fit.log path using known structures.

    Recognised structures:
      - EOS output: <eos_output>/<ConvinoConfigName>/<CorrelationName>/<Value>/.../fit.log
      - Job folder: <jobs_folder>/<CorrelationName>/<ConvinoConfigName__CorrelationName__Value>/fit.log

    Returns (scanname, value) where value is None for nominal or unknown.
    """
    parts = [p for p in path.split(os.sep) if p != '']
    scanname = 'unknown'
    value = None

    # check EOS structure first (more deterministic)
    if eos_output:
        eos_parts = [p for p in os.path.abspath(eos_output).split(os.sep) if p != '']
        idx = _find_sublist_index(parts, eos_parts)
        if idx >= 0:
            rel = parts[idx + len(eos_parts):]
            # expected rel: [ConvinoConfigName?, CorrelationName, Value, ...]
            if len(rel) > 0:
                # try to find the last numeric token in the relative tail
                val = None
                for i in range(len(rel) - 1, -1, -1):
                    try:
                        val = float(rel[i])
                        val_idx = i
                        break
                    except Exception:
                        val = None
                if val is not None:
                    value = val
                    # prefer the token before the value as scanname if it looks like a correlation name
                    if val_idx - 1 >= 0:
                        scanname = rel[val_idx - 1]
                    elif len(rel) >= 2:
                        scanname = rel[0]
                    return scanname, value
                else:
                    # no numeric token found: try to pick correlation name as second element if present
                    if len(rel) >= 2:
                        scanname = rel[1]
                    else:
                        scanname = rel[0]
                    return scanname, None

    # check job-folder structure
    if jobs_folder:
        jb_parts = [p for p in os.path.abspath(jobs_folder).split(os.sep) if p != '']
        idx = _find_sublist_index(parts, jb_parts)
        if idx >= 0:
            rel = parts[idx + len(jb_parts):]
            # expect rel: [CorrelationName, JobFolder, ...]
            if len(rel) >= 1:
                scanname = rel[0]
            if len(rel) >= 2:
                jobfolder = rel[1]
                if jobfolder.lower() == 'nominal':
                    return scanname, None
                # try to extract value after final '__' in jobfolder
                if '__' in jobfolder:
                    try:
                        valstr = jobfolder.rsplit('__', 1)[1]
                        value = float(valstr)
                        return scanname, value
                    except Exception:
                        value = None
            # if not found in jobfolder, search rel for last numeric token
            for i in range(len(rel) - 1, -1, -1):
                try:
                    v = float(rel[i])
                    return scanname, v
                except Exception:
                    continue

    # fallback: look for a folder that looks like 'lhs__rhs' for scanname
    for p in reversed(parts):
        if '__' in p and not p.endswith('.tgz'):
            scanname = p
            break

    # fallback: find any numeric token as value
    for p in reversed(parts):
        try:
            v = float(p)
            value = v
            break
        except Exception:
            if '__' in p:
                try:
                    tok = p.rsplit('__', 1)[1]
                    value = float(tok)
                    break
                except Exception:
                    pass

    return scanname, value


def collect_results(jobs_folder: str, eos_output: Optional[str], verbose: bool = False, debug: bool = False) -> Dict[str, dict]:
    """Collect all parsed fit.log results and organize by scan.

    Returned structure:
      { scanname: { 'entries': [ {value, central, exp_unc, pdf_unc, total_unc, path}, ... ],
                     'nominal': { ... } or None } }
    """
    # bases = [jobs_folder]
    bases = []
    if eos_output:
        bases.append(eos_output)

    # find all convino_run.log files under the bases; convino_run.log always exists per-job using find_convino_logs
    convino_logs = find_convino_logs(bases)

    if verbose:
        print(f'Found {len(convino_logs)} convino_run.log files under: {bases}')
    if debug:
        for p in convino_logs:
            print('  -', p)

    scans = defaultdict(lambda: {'entries': [], 'nominal': None})

    for convino in convino_logs:
        # corresponding fit.log lives in same directory if invertible; check it first
        convino_path = convino
        f = os.path.join(os.path.dirname(convino_path), 'fit.log')
        parsed = None
        noninvert = False
        # if fit.log exists, parse it and assume invertible (skip reading convino.log)
        if os.path.isfile(f):
            parsed = parse_fit_log(f)
        else:
            # no fit.log -> inspect convino.log for non-invertible marker
            try:
                txt = open(convino_path, 'r', errors='ignore').read()
                if 'RobustInvert: Cholesky failed' in txt:
                    noninvert = True
            except Exception:
                if debug:
                    print('DEBUG: failed to read', convino_path)

        # if there is neither a parsed fit nor a non-invertible flag, skip
        if parsed is None and not noninvert:
            if debug:
                print('DEBUG: no fit result and not non-invertible for', convino_path)
            continue

        scanname, value = infer_scan_and_value_from_path(f if os.path.isfile(f) else convino_path, jobs_folder, eos_output)

        entry = {
            'value': value,
            'central': (parsed['central'] if parsed is not None else None),
            'exp_unc': (parsed['exp_unc'] if parsed is not None else None),
            'pdf_unc': (parsed['pdf_unc'] if parsed is not None else None),
            'total_unc': (parsed['total_unc'] if parsed is not None else None),
            'path': (parsed['path'] if parsed is not None else os.path.abspath(convino_path)),
            'nonInvertible': noninvert,
            'convino_log': convino_path,
        }
        if verbose or debug:
            print(f"Parsed/Found: convino={convino_path}, fit={f if os.path.isfile(f) else 'N/A'}\n  -> scan={scanname}, value={value}, central={entry['central']}, total_unc={entry['total_unc']}, nonInvertible={noninvert}")

        key_path = f if os.path.isfile(f) else convino_path
        if (os.path.sep + 'nominal' + os.path.sep) in key_path or os.path.basename(os.path.dirname(key_path)).lower().startswith('nominal'):
            if verbose:
                print(f"Marking nominal for scan {scanname} from {key_path}")
            scans[scanname]['nominal'] = entry
        else:
            scans[scanname]['entries'].append(entry)

    # sort entries by value
    for s in scans:
        scans[s]['entries'] = sorted(scans[s]['entries'], key=lambda e: (float('inf') if e['value'] is None else e['value']))

    return scans


def plot_scan(scanname: str, scan: dict, outdir: str, unblind: bool = False) -> None:
    """Create a 3-panel plot for a single scan and save PNG."""
    entries = scan.get('entries', [])
    nominal = scan.get('nominal')
    if len(entries) == 0:
        return

    xs = np.array([e['value'] if e['value'] is not None else np.nan for e in entries], dtype=float)
    centrals = np.array([e['central'] for e in entries], dtype=float)
    total_uncs = np.array([e['total_unc'] for e in entries], dtype=float)

    # since there are some weird lines in the plot, pls sort them by the xs
    sort_idx = np.argsort(xs)
    xs = xs[sort_idx]
    centrals = centrals[sort_idx]
    total_uncs = total_uncs[sort_idx]
    # reorder entries to match sorted arrays so we can color points consistently
    try:
        sort_list = list(map(int, sort_idx.tolist())) if hasattr(sort_idx, 'tolist') else list(map(int, sort_idx))
    except Exception:
        sort_list = list(sort_idx)
    entries_sorted = [entries[i] for i in sort_list]

    if nominal is not None:
        nom_c = nominal['central']
        nom_unc = nominal['total_unc']
    else:
        # fall back to first entry as nominal if absent
        nom_c = centrals[0]
        nom_unc = total_uncs[0]

    # mask of entries that are exactly identical to nominal central
    try:
        same_as_nom = np.array([(e.get('central') is not None and nom_c is not None and e.get('central') == nom_c) for e in entries_sorted], dtype=bool)
    except Exception:
        same_as_nom = np.zeros(len(entries_sorted), dtype=bool)
    # check if there are multiple matches. If so, and the name includes "SinglePart" choose the point where the correlation value is 0.75
    # if 'SinglePart' in scanname:
        # print ('DEBUG: same_as_nom before adjustment for SinglePart:', same_as_nom)
    if np.sum(same_as_nom) > 1 and 'SinglePart' in scanname:
        for i, e in enumerate(entries_sorted):
            if e.get('value') == 0.75:
                same_as_nom = np.zeros(len(entries_sorted), dtype=bool)
                same_as_nom[i] = True
                break
            # ignore all values for the central that are not 0.0, abs(0.25), abs(0.5), abs(0.75)
            if e.get('value') not in [0.0, 0.25, 0.5, 0.75, -0.25, -0.5, -0.75]:
                same_as_nom[i] = False

    # color array: red for exactly-nominal entries, default blue otherwise
    colors = np.where(same_as_nom, 'red', 'C0')

    delta_c = centrals - nom_c
    # relative deviation wrt nominal central value (fractional)
    try:
        if nom_c != 0:
            rel_dev = delta_c / nom_c
        else:
            rel_dev = np.full_like(delta_c, np.nan)
    except Exception:
        rel_dev = np.full_like(delta_c, np.nan)
    rel_unc = total_uncs / centrals
    rel_unc_vs_nom = total_uncs / nom_unc

    # nominal relative (nom_unc / nom_c) precomputed for plotting
    try:
        nom_rel = nom_unc / nom_c if (nom_c and not math.isnan(nom_c)) else float('nan')
    except Exception:
        nom_rel = float('nan')

    # try to apply mplhep CMS style if available
    try:
        import mplhep as hep

        plt.style.use(hep.style.CMS)
    except Exception:
        pass

    # layout: 3 rows x 2 cols GridSpec
    from matplotlib import gridspec

    fig = plt.figure(figsize=(10, 10))
    # one large title on top with the scanname
    fig.suptitle(f'Scan summary: {scanname}', fontsize=12)
    # reduce horizontal spacing between columns so panels sit closer
    gs = gridspec.GridSpec(3, 2, width_ratios=[1, 1], height_ratios=[1, 1, 1], hspace=0.25, wspace=0.18)

    # Left column: two panels (central vs correlation, total_unc vs correlation)
    ax_l_top = fig.add_subplot(gs[0, 0])
    ax_l_mid = fig.add_subplot(gs[1, 0])
    # Bottom-left: relative deviation (Fitted - nominal) / nominal
    ax_l_bot = fig.add_subplot(gs[2, 0])

    # Right column: three panels (delta central, relative uncertainty, uncertainty/nominal)
    ax_r_top = fig.add_subplot(gs[0, 1])
    ax_r_mid = fig.add_subplot(gs[1, 1])
    ax_r_bot = fig.add_subplot(gs[2, 1])

    # reduce font sizes for labels and ticks
    fs_label = 8
    fs_ticks = 7

    # define masks requiring both x and y to be finite for each plot to avoid NaN limits
    mask_x = ~np.isnan(xs)
    mask_c = ~np.isnan(centrals)
    mask_tu = ~np.isnan(total_uncs)
    mask_rel_dev = ~np.isnan(rel_dev)
    mask_rel_unc = ~np.isnan(rel_unc)
    mask_rel_unc_vs_nom = ~np.isnan(rel_unc_vs_nom)

    mask_top = mask_x & mask_c
    mask_mid = mask_x & mask_tu
    mask_bot = mask_x & mask_rel_dev

    # top-left: central vs correlation
    if np.any(mask_top):
        if unblind:
            ax_l_top.scatter(xs[mask_top], centrals[mask_top], c=colors[mask_top], marker='o', s=30, label='Scan points', edgecolors='k')
            cmin = np.nanmin(centrals[mask_top])
            cmax = np.nanmax(centrals[mask_top])
            desired_top_l_min = nom_c - 0.05
            desired_top_l_max = nom_c + 0.05
            if (cmin >= desired_top_l_min) and (cmax <= desired_top_l_max):
                ax_l_top.set_ylim(desired_top_l_min, desired_top_l_max)
            else:
                crange = cmax - cmin
                pad = 0.02 * crange if crange > 0 else max(1e-3, 0.01 * abs(cmin))
                ax_l_top.set_ylim(cmin - pad, cmax + pad)
            # format y axis to plain numbers (no scientific notation)
            from matplotlib.ticker import ScalarFormatter

            sf = ScalarFormatter()
            sf.set_scientific(False)
            sf.set_useOffset(False)
            ax_l_top.yaxis.set_major_formatter(sf)
            # draw nominal horizontal line for fitted mass
            try:
                ax_l_top.axhline(nom_c, color='red', linestyle='--', lw=0.5, label='Nominal')
            except Exception:
                pass
        else:
            ax_l_top.text(0.0, 0.5, 'Blinded.', ha='center', va='center', fontsize=fs_label)
    else:
        ax_l_top.text(0.5, 0.5, 'no correlation values', ha='center', va='center', fontsize=fs_label)

    # middle-left: total uncertainty vs correlation
    if np.any(mask_mid):
        ax_l_mid.scatter(xs[mask_mid], total_uncs[mask_mid], c=colors[mask_mid], marker='o', s=30, label='Scan points', edgecolors='k')
        try:
            ax_l_mid.axhline(nom_unc, color='red', linestyle='--', lw=0.5, label='Nominal')
        except Exception:
            pass
    else:
        ax_l_mid.text(0.5, 0.5, 'no correlation values', ha='center', va='center', fontsize=fs_label)

    # bottom-left: relative deviation
    if np.any(mask_bot):
        ax_l_bot.scatter(xs[mask_bot], rel_dev[mask_bot], c=colors[mask_bot], marker='o', s=30, label='Scan points', edgecolors='k')
        try:
            ax_l_bot.axhline(0.0, color='red', linestyle='--', lw=0.5, label='Nominal')
        except Exception:
            pass
    else:
        ax_l_bot.text(0.5, 0.5, 'no correlation values', ha='center', va='center', fontsize=fs_label)

    # right-middle: relative uncertainty (plot even if some values missing)
    if np.any(mask_x & mask_rel_unc):
        ax_r_mid.scatter(xs[mask_x & mask_rel_unc], rel_unc[mask_x & mask_rel_unc], c=colors[mask_x & mask_rel_unc], marker='o', s=30, label='Scan points', edgecolors='k')
    else:
        ax_l_top.text(0.5, 0.5, 'no correlation values', ha='center', va='center', fontsize=fs_label)
        ax_l_mid.text(0.5, 0.5, 'no correlation values', ha='center', va='center', fontsize=fs_label)


    ax_l_top.set_ylabel(r'Fitted $m_{t}$ [GeV]', fontsize=fs_label)
    ax_l_top.set_xlabel('Correlation', fontsize=fs_label)
    ax_l_mid.set_ylabel('Total uncertainty [GeV]', fontsize=fs_label)
    ax_l_mid.set_xlabel('Correlation', fontsize=fs_label)
    ax_l_bot.set_ylabel('Relative deviation', fontsize=fs_label)
    ax_l_bot.set_xlabel('Correlation', fontsize=fs_label)
    ax_l_top.tick_params(axis='both', labelsize=fs_ticks)
    ax_l_mid.tick_params(axis='both', labelsize=fs_ticks)
    ax_l_bot.tick_params(axis='both', labelsize=fs_ticks)

    # set range always from -1.1 to 1.1 for correlation x-axis
    ax_l_top.set_xlim(-1.1, 1.1)
    ax_l_mid.set_xlim(-1.1, 1.1)
    ax_l_bot.set_xlim(-1.1, 1.1)

    # detect entries that are non-invertible and draw dashed gray vertical areas
    try:
        noninv_flags = np.array([bool(e.get('nonInvertible')) for e in entries], dtype=bool)
        noninv_flags = noninv_flags[sort_idx]
        # width for shaded area around the x position
        w = 0.02
        for xi, ni in zip(xs, noninv_flags):
            if not ni or np.isnan(xi):
                continue
            for _ax in (ax_l_mid, ax_l_bot, ax_r_top, ax_r_mid, ax_r_bot):
                try:
                    _ax.axvspan(xi - w, xi + w, facecolor='gray', alpha=0.25, edgecolor='gray', linestyle='--', linewidth=0.8, zorder=0)
                except Exception:
                    pass
            if unblind:
                try:
                    ax_l_top.axvspan(xi - w, xi + w, facecolor='gray', alpha=0.25, edgecolor='gray', linestyle='--', linewidth=0.8, zorder=0)
                except Exception:
                    pass
    except Exception:
        pass

    # draw linear fits (degree=1) on all six subplots if there are >=2 valid points
    x_line = np.linspace(-1, 1, 200)
    def _fit_and_plot(ax, xvals, yvals):
        try:
            mask2 = (~np.isnan(xvals)) & (~np.isnan(yvals))
            if np.sum(mask2) >= 2:
                coeff = np.polyfit(xvals[mask2], yvals[mask2], 2)
                y_line = np.polyval(coeff, x_line)
                ax.plot(x_line, y_line, linestyle='--', color='blue', linewidth=1, alpha=0.7, label='Fit')
        except Exception:
            pass

    if unblind:
        _fit_and_plot(ax_l_top, xs, centrals)
    _fit_and_plot(ax_l_mid, xs, total_uncs)
    _fit_and_plot(ax_l_bot, xs, rel_dev)
    _fit_and_plot(ax_r_top, xs, delta_c)
    _fit_and_plot(ax_r_mid, xs, rel_unc)
    _fit_and_plot(ax_r_bot, xs, rel_unc_vs_nom)

    # Y-axis zoom constraints (allow zooming out if data exceeds these ranges)
    # middle-left (total uncertainty): nominal +/- 0.05
    try:
        tu_min = np.nanmin(total_uncs)
        tu_max = np.nanmax(total_uncs)
        desired_mid_l_min = nom_unc - 0.05
        desired_mid_l_max = nom_unc + 0.05
        if (tu_min >= desired_mid_l_min) and (tu_max <= desired_mid_l_max):
            ax_l_mid.set_ylim(desired_mid_l_min, desired_mid_l_max)
        else:
            trange = tu_max - tu_min
            tpad = 0.02 * trange if trange > 0 else max(1e-6, 0.01 * abs(tu_min))
            ax_l_mid.set_ylim(tu_min - tpad, tu_max + tpad)
    except Exception:
        pass

    # bottom-left (relative deviation): fixed range [-0.01, 0.01], but expand if data outside
    try:
        rd_min = np.nanmin(rel_dev)
        rd_max = np.nanmax(rel_dev)
        if (rd_min >= -0.01) and (rd_max <= 0.01):
            ax_l_bot.set_ylim(-0.01, 0.01)
        else:
            rrange = rd_max - rd_min
            rpad = 0.02 * rrange if rrange > 0 else 1e-6
            ax_l_bot.set_ylim(rd_min - rpad, rd_max + rpad)
    except Exception:
        pass

    # right column plots
    # right-top: color points red when exactly nominal
    mask_dc = (~np.isnan(xs)) & (~np.isnan(delta_c))
    if np.any(mask_dc):
        ax_r_top.scatter(xs[mask_dc], delta_c[mask_dc], c=colors[mask_dc], marker='o', s=30, label='Scan points', edgecolors='k')
    ax_r_top.axhline(0, color='red', linestyle='--', lw=0.5, label='Nominal')
    ax_r_top.set_ylabel('Fitted - nominal $m_{t}$ [GeV]', fontsize=fs_label)
    ax_r_top.set_xlabel('Correlation', fontsize=fs_label)
    ax_r_top.tick_params(axis='both', labelsize=fs_ticks)

    # ax_r_mid.plot(xs, rel_unc, marker='o')
    ax_r_mid.set_ylabel('Relative uncertainty', fontsize=fs_label)
    ax_r_mid.set_xlabel('Correlation', fontsize=fs_label)
    ax_r_mid.tick_params(axis='both', labelsize=fs_ticks)
    # draw nominal horizontal line for relative uncertainty (nom_unc/nom_c)
    try:
        if not math.isnan(nom_rel):
            ax_r_mid.axhline(nom_rel, color='red', linestyle='--', lw=0.5, label='Nominal')
    except Exception:
        pass

    mask_rbot = (~np.isnan(xs)) & (~np.isnan(rel_unc_vs_nom))
    if np.any(mask_rbot):
        ax_r_bot.scatter(xs[mask_rbot], rel_unc_vs_nom[mask_rbot], c=colors[mask_rbot], marker='o', s=30, label='Scan points', edgecolors='k')
    ax_r_bot.axhline(1., color='red', linestyle='--', lw=0.5, label='Nominal')
    ax_r_bot.set_ylabel('Uncertainty / nominal uncertainty', fontsize=fs_label)
    ax_r_bot.set_xlabel('Correlation', fontsize=fs_label)
    ax_r_bot.tick_params(axis='both', labelsize=fs_ticks)

    # set range always from -1.1 to 1.1 for correlation x-axis
    ax_r_top.set_xlim(-1.1, 1.1)
    ax_r_mid.set_xlim(-1.1, 1.1)
    ax_r_bot.set_xlim(-1.1, 1.1)

    # right column y-limits with allowed expansion
    # top-right (delta central): prefer [-0.05, 0.05]
    try:
        dc_min = np.nanmin(delta_c)
        dc_max = np.nanmax(delta_c)
        if (dc_min >= -0.05) and (dc_max <= 0.05):
            ax_r_top.set_ylim(-0.05, 0.05)
        else:
            drange = dc_max - dc_min
            dpad = 0.02 * drange if drange > 0 else 1e-6
            ax_r_top.set_ylim(dc_min - dpad, dc_max + dpad)
    except Exception:
        pass

    # middle-right (relative uncertainty): nominal +/- 0.0001 (nom_rel = nom_unc/nom_c)
    try:
        nom_rel = nom_unc / nom_c if (nom_c and not math.isnan(nom_c)) else float('nan')
        ru_min = np.nanmin(rel_unc)
        ru_max = np.nanmax(rel_unc)
        desired_mid_r_min = nom_rel - 0.0001
        desired_mid_r_max = nom_rel + 0.0001
        if (ru_min >= desired_mid_r_min) and (ru_max <= desired_mid_r_max):
            ax_r_mid.set_ylim(desired_mid_r_min, desired_mid_r_max)
        else:
            rrange = ru_max - ru_min
            rpad = 0.02 * rrange if rrange > 0 else max(1e-6, 0.01 * abs(ru_min))
            ax_r_mid.set_ylim(ru_min - rpad, ru_max + rpad)
    except Exception:
        pass

    # bottom-right (uncertainty / nominal uncertainty): prefer [0.999, 1.001]
    try:
        rv_min = np.nanmin(rel_unc_vs_nom)
        rv_max = np.nanmax(rel_unc_vs_nom)
        if (rv_min >= 0.99) and (rv_max <= 1.01):
            ax_r_bot.set_ylim(0.99, 1.01)
        else:
            vrange = rv_max - rv_min
            vpad = 0.02 * vrange if vrange > 0 else 1e-6
            ax_r_bot.set_ylim(rv_min - vpad, rv_max + vpad)
    except Exception:
        pass

    # add deduplicated legends to each axis
    def _add_legend(ax):
        try:
            handles, labels = ax.get_legend_handles_labels()
            if not handles:
                return
            seen = set()
            new_h = []
            new_l = []
            for h, l in zip(handles, labels):
                if l in seen:
                    continue
                seen.add(l)
                new_h.append(h)
                new_l.append(l)
            if new_h:
                ax.legend(new_h, new_l, fontsize=8, loc='best', frameon=False)
        except Exception:
            pass

    for _ax in (ax_l_mid, ax_l_bot, ax_r_top, ax_r_mid, ax_r_bot):
        _add_legend(_ax)
    if unblind:
        _add_legend(ax_l_top)

    os.makedirs(outdir, exist_ok=True)
    outfn_png = os.path.join(outdir, f'scan_{scanname}.png')
    outfn_pdf = os.path.join(outdir, f'scan_{scanname}.pdf')
    # reduce outer margins to give more room for the two columns
    fig.subplots_adjust(top=0.95, bottom=0.05, left=0.06, right=0.97)
    fig.savefig(outfn_png)
    try:
        fig.savefig(outfn_pdf)
    except Exception:
        pass
    plt.close(fig)


def normalize_scans(scans: Dict[str, dict], verbose: bool = False) -> None:
    """Ensure each scan has a `nominal` entry.

    Rules:
      - If any entry in `entries` has a path whose basename is 'fit.log' or
        has value==None, and `nominal` is empty, promote that entry to `nominal`.
      - If `nominal` exists but is not present in `entries`, append it to `entries`.
    This modifies `scans` in-place.
    """
    for sname, sdict in scans.items():
        if not isinstance(sdict, dict):
            continue
        entries = sdict.get('entries', []) or []
        nominal = sdict.get('nominal')

        # find candidate in entries
        if nominal is None:
            cand_idx = None
            for i, e in enumerate(entries):
                p = e.get('path', '') if isinstance(e, dict) else ''
                val = e.get('value') if isinstance(e, dict) else None
                if os.path.basename(p) == 'fit.log' or val is None:
                    cand_idx = i
                    break
            if cand_idx is not None:
                sdict['nominal'] = entries.pop(cand_idx)
                if verbose:
                    print(f"normalize_scans: promoted nominal for {sname} from entries[{cand_idx}]")

        # ensure nominal is present in entries as well
        nominal = sdict.get('nominal')
        if nominal is not None:
            paths = [e.get('path', '') for e in entries if isinstance(e, dict)]
            npath = nominal.get('path', '') if isinstance(nominal, dict) else ''
            if npath and (npath not in paths):
                entries.append(nominal)
                sdict['entries'] = entries
                if verbose:
                    print(f"normalize_scans: appended nominal to entries for {sname}")


def propagate_global_nominal(scans: Dict[str, dict], verbose: bool = False) -> None:
    """If a top-level scan named 'nominal' exists, use its nominal entry as the
    common nominal and set it for all other scans that don't already have one.
    """
    global_scan = scans.get('nominal') or scans.get('Nominal') or scans.get('fit.log')
    if not global_scan:
        if verbose:
            print("propagate_global_nominal: no global nominal scan found")
        return

    # find the representative nominal entry in the global_scan
    rep = None
    if isinstance(global_scan, dict):
        rep = global_scan.get('nominal') or (global_scan.get('entries')[0] if global_scan.get('entries') else None)

    if not rep:
        if verbose:
            print("propagate_global_nominal: no nominal entry found in global nominal scan")
        return

    for sname, sdict in scans.items():
        if sname == 'nominal':
            continue
        if not isinstance(sdict, dict):
            continue
        cur_nom = sdict.get('nominal')
        # if cur_nom is None:
        sdict['nominal'] = rep
        # also ensure present in entries
        entries = sdict.get('entries') or []
        paths = [e.get('path', '') for e in entries if isinstance(e, dict)]
        if rep.get('path', '') not in paths:
            entries.append(rep)
            sdict['entries'] = entries
        if verbose:
            print(f"propagate_global_nominal: set nominal for {sname} from global nominal")


def plot_2d(scans: Dict[str, dict], outdir: str):
    """Create two 2D summary plots across scans.

    - Z1: central difference vs nominal
    - Z2: total_unc - nominal_total_unc
    """
    names = sorted(scans.keys())
    # fill the names with dummy values from 0 to N to make shorter labels
    namesNew = [str(i) for i in range(len(names))]

    # collect unique sorted values across all scans
    all_vals = sorted({v for s in scans.values() for e in s.get('entries', []) for v in ([e['value']] if e['value'] is not None else [])})
    if len(all_vals) == 0:
        return

    val_to_idx = {v: i for i, v in enumerate(all_vals)}

    Zc = np.full((len(all_vals), len(names)), np.nan)
    Zu = np.full((len(all_vals), len(names)), np.nan)
    NonInv = np.zeros((len(all_vals), len(names)), dtype=bool)
    Same = np.zeros((len(all_vals), len(names)), dtype=bool)

    for j, name in enumerate(names):
        s = scans[name]
        nom = s.get('nominal')
        nom_c = nom['central'] if nom else np.nan
        nom_u = nom['total_unc'] if nom else np.nan

        for e in s.get('entries', []):
            v = e['value']
            if v is None:
                continue
            i = val_to_idx.get(v)
            if i is None:
                continue
            # print (j, name, v, i, e['central'], nom_c)
            # for non-invertible entries, the values may be None
            if e['central'] is None or e['total_unc'] is None:
                if e.get('nonInvertible'):
                    NonInv[i, j] = True
                continue
            Zc[i, j] = e['central'] - nom_c
            Zu[i, j] = e['total_unc'] - nom_u
            if e.get('nonInvertible'):
                NonInv[i, j] = True
            # mark if this entry's central equals the nominal central exactly
            try:
                if (e.get('central') is not None) and (nom_c is not None) and (e.get('central') == nom_c):
                    Same[i, j] = True
            except Exception:
                pass
        # check if there are multiple matches. If so, and the name includes "SinglePart" choose the point where the correlation value is 0.75
        # if 'SinglePart' in name:
        if True:
            # print ('DEBUG: Same[:, j] before adjustment for SinglePart:', Same[:, j])
            count_same = np.sum(Same[:, j])
            if count_same > 1:
                for i_val, e in enumerate(s.get('entries', [])):
                    if e.get('value') == 0.75:
                        Same[:, j] = False
                        idx = val_to_idx.get(e.get('value'))
                        if idx is not None:
                            Same[idx, j] = True
                        break
                    # if any value not in [0.0, abs(0.25), abs(0.5), abs(0.75)], set Same to False
                    if e.get('value') not in [0.0, 0.25, 0.5, 0.75, -0.25, -0.5, -0.75]:
                        idx = val_to_idx.get(e.get('value'))
                        if idx is not None:
                            Same[idx, j] = False

    os.makedirs(outdir, exist_ok=True)

    # plot central differences
    fig, ax = plt.subplots(1, 1, figsize=(max(6, len(namesNew) * 0.4), max(4, len(all_vals) * 0.4)))
    # center the z axis around 0
    im = ax.imshow(Zc, aspect='auto', origin='lower', interpolation='none', vmin=-np.nanmax(np.abs(Zc)), vmax=np.nanmax(np.abs(Zc)))
    ax.set_xticks(range(len(namesNew)))
    ax.set_xticklabels(namesNew, rotation=90)
    ax.set_yticks(range(len(all_vals)))
    ax.set_yticklabels([f'{v:.3f}' for v in all_vals])
    ax.set_xlabel('Scan index')
    ax.set_ylabel('Correlation')
    # ax.set_title(r'Fitted - nominal $m_{t}$ [GeV]')
    # instead of a title, put a label on the z axis
    fig.colorbar(im, ax=ax, label=r'Fitted - nominal $m_{t}$ [GeV]')
    # avoid tight_layout warning for mixed Axes; use subplots_adjust and expand width
    fig.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.12)
    # overlay dashed gray boxes for non-invertible points
    try:
        from matplotlib.patches import Rectangle

        for (i, j), flag in np.ndenumerate(NonInv):
            if not flag:
                continue
            # rectangle centered at (j, i) with size 1x1 in image coordinates
            rect = Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, facecolor='gray', edgecolor='gray', linestyle='-', linewidth=1, fill = True, hatch='//', alpha=0.4)
            ax.add_patch(rect)
        # overlay red outline for entries exactly equal to nominal
        for (i, j), flag in np.ndenumerate(Same):
            if not flag:
                continue
            try:
                rect = Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, facecolor='none', edgecolor='red', linestyle='-', linewidth=1.4)
                ax.add_patch(rect)
            except Exception:
                pass
    except Exception:
        pass
    # put a small label outside of the plot explaining what the hatched area means
    # ax.text(1.02, 0.02, 'Hatched: non-invertible', transform=ax.transAxes, fontsize=8, va='bottom', ha='left', color='gray')
    hep.cms.label(exp="ATLAS+CMS", llabel= "Work in Progress", rlabel = "8+13 TeV", ax=ax)
    p_png = os.path.join(outdir, '2d_central_diff.png')
    p_pdf = os.path.join(outdir, '2d_central_diff.pdf')
    fig.savefig(p_png)
    try:
        fig.savefig(p_pdf)
    except Exception:
        pass
    plt.close(fig)

    # and the second plot
    fig, ax = plt.subplots(1, 1, figsize=(max(6, len(namesNew) * 0.4), max(4, len(all_vals) * 0.4)))
    # center the z axis around 0
    im = ax.imshow(Zu, aspect='auto', origin='lower', interpolation='none', vmin=-np.nanmax(np.abs(Zu)), vmax=np.nanmax(np.abs(Zu)))
    ax.set_xticks(range(len(namesNew)))
    ax.set_xticklabels(namesNew, rotation=90)
    ax.set_yticks(range(len(all_vals)))
    ax.set_yticklabels([f'{v:.3f}' for v in all_vals])
    ax.set_xlabel('Scan index')
    ax.set_ylabel('Correlation')
    # ax.set_title('Uncertainty - nominal uncertainty [GeV]')
    fig.colorbar(im, ax=ax, label='Uncertainty - nominal uncertainty [GeV]')
    # avoid tight_layout warning for mixed Axes; use subplots_adjust and expand width
    fig.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.12)
    # overlay dashed gray boxes for non-invertible points
    try:
        from matplotlib.patches import Rectangle

        for (i, j), flag in np.ndenumerate(NonInv):
            if not flag:
                continue
            # rectangle centered at (j, i) with size 1x1 in image coordinates
            rect = Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, facecolor='gray', edgecolor='gray', linestyle='-', linewidth=1, fill = True, hatch='//', alpha=0.4)
            ax.add_patch(rect)
        # overlay red outline for entries exactly equal to nominal
        for (i, j), flag in np.ndenumerate(Same):
            if not flag:
                continue
            try:
                rect = Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, facecolor='none', edgecolor='red', linestyle='-', linewidth=1.4)
                ax.add_patch(rect)
            except Exception:
                pass
    except Exception:
        pass
    # put a small label outside of the plot explaining what the hatched area means
    # ax.text(1.02, 0.02, 'Hatched: non-invertible', transform=ax.transAxes, fontsize=8, va='bottom', ha='left', color='gray')
    hep.cms.label(exp="ATLAS+CMS", llabel= "Work in Progress", rlabel = "8+13 TeV", ax=ax)
    p_png = os.path.join(outdir, '2d_total_unc_diff.png')
    p_pdf = os.path.join(outdir, '2d_total_unc_diff.pdf')
    fig.savefig(p_png)
    try:
        fig.savefig(p_pdf)
    except Exception:
        pass
    plt.close(fig)


def compute_scan_correlations(scans: Dict[str, dict], verbose: bool = False) -> List[dict]:
    """Compute linear correlation (Pearson r) between mass central value and correlation for each scan.

    Returns list of dicts with keys: scan, r, slope, intercept, n_points
    sorted by absolute r desc.
    """
    rows = []
    for name, s in scans.items():
        if name.lower() == 'nominal':
            continue
        if not isinstance(s, dict):
            continue
        entries = s.get('entries', []) or []
        xs = []
        ys = []
        for e in entries:
            v = e.get('value')
            c = e.get('central')
            try:
                if v is None or c is None:
                    continue
                xv = float(v)
                yv = float(c)
                xs.append(xv)
                ys.append(yv)
            except Exception:
                continue
        n = len(xs)
        if verbose:
            print(f'DEBUG: scan {name} has {n} valid points for correlation computation')
            print('  xs =', xs)
            print('  ys =', ys)
        if n < 2:
            continue
        xs_arr = np.array(xs)
        ys_arr = np.array(ys)
        # Pearson r
        try:
            r = float(np.corrcoef(xs_arr, ys_arr)[0, 1])
        except Exception:
            r = float('nan')
        # linear fit
        try:
            slope, intercept = np.polyfit(xs_arr, ys_arr, 2)
        except Exception:
            slope, intercept = float('nan'), float('nan')

        rows.append({'scan': name, 'r': r, 'slope': slope, 'intercept': intercept, 'n_points': n})

    rows = sorted(rows, key=lambda r: abs(r['slope']) if not math.isnan(r['slope']) else 0.0, reverse=True)
    return rows


def print_top_deviations(scans: Dict[str, dict], topn: int = 10, unblind: bool = False) -> None:
    """Print top-N scans ranked by maximum absolute central deviation w.r.t. nominal.

    For each scan, compute max |central - nominal_c| across entries. Report
    the correlation value where it occurs and the ratio to the nominal uncertainty.
    """
    rows = []
    for name, s in scans.items():
        if name.lower() == 'nominal':
            continue
        if not isinstance(s, dict):
            continue
        nominal = s.get('nominal')
        entries = s.get('entries', []) or []
        if not nominal or len(entries) == 0:
            continue
        nom_c = nominal.get('central')
        nom_u = nominal.get('total_unc') if nominal.get('total_unc') is not None else float('nan')

        max_delta = 0.0
        best = None
        for e in entries:
            try:
                delta = abs(e.get('central', 0.0) - nom_c)
            except Exception:
                continue
            if delta > max_delta:
                max_delta = delta
                best = e

        if best is None:
            continue

        ratio = (max_delta / nom_u) if (nom_u and not math.isnan(nom_u) and nom_u != 0) else None
        rows.append({'scan': name, 'max_abs_delta': max_delta, 'value': best.get('value') if unblind else max_delta, 'path': best.get('path'), 'ratio_nominal_unc': ratio})

    rows = sorted(rows, key=lambda r: r['max_abs_delta'], reverse=True)
    print('\nTop {} scans by max absolute central deviation (w.r.t. nominal):'.format(min(topn, len(rows))))
    print('{:3s} {:60s} {:>10s} {:>10s} {:>10s}'.format('#', 'scan', 'maxDelta[GeV]', 'corr', 'ratio'))
    for i, r in enumerate(rows[:topn]):
        ratio_str = f"{r['ratio_nominal_unc']:.3f}" if r['ratio_nominal_unc'] is not None else 'n/a'
        val_str = f"{r['value']:.3f}" if (r['value'] is not None) else 'n/a'
        print(f"{i+1:3d} {r['scan'][:60]:60s} {r['max_abs_delta']:10.4f} {val_str:>10s} {ratio_str:>10s}")
    return rows


def print_top_uncertainty_deviations(scans: Dict[str, dict], topn: int = 10) -> None:
    """Print top-N scans ranked by maximum absolute uncertainty deviation (total_unc - nominal_total_unc).

    For each scan, compute max |total_unc - nominal_total_unc| across entries. Report
    the correlation value where it occurs and the ratio to the nominal uncertainty.
    """
    rows = []
    for name, s in scans.items():
        if name.lower() == 'nominal':
            continue
        if not isinstance(s, dict):
            continue
        nominal = s.get('nominal')
        entries = s.get('entries', []) or []
        if not nominal or len(entries) == 0:
            continue
        nom_u = nominal.get('total_unc') if nominal.get('total_unc') is not None else float('nan')

        max_du = 0.0
        best = None
        for e in entries:
            try:
                du = abs(e.get('total_unc', 0.0) - nom_u)
            except Exception:
                continue
            if du > max_du:
                max_du = du
                best = e

        if best is None:
            continue

        ratio = (max_du / nom_u) if (nom_u and not math.isnan(nom_u) and nom_u != 0) else None
        rows.append({'scan': name, 'max_abs_unc_diff': max_du, 'value': best.get('value'), 'path': best.get('path'), 'ratio_nominal_unc': ratio})

    rows = sorted(rows, key=lambda r: r['max_abs_unc_diff'], reverse=True)
    print('\nTop {} scans by max absolute uncertainty deviation (w.r.t. nominal):'.format(min(topn, len(rows))))
    print('{:3s} {:60s} {:>12s} {:>10s} {:>10s}'.format('#', 'scan', 'maxDeltaUnc[GeV]', 'corr', 'ratio'))
    for i, r in enumerate(rows[:topn]):
        ratio_str = f"{r['ratio_nominal_unc']:.3f}" if r['ratio_nominal_unc'] is not None else 'n/a'
        val_str = f"{r['value']:.3f}" if (r['value'] is not None) else 'n/a'
        print(f"{i+1:3d} {r['scan'][:60]:60s} {r['max_abs_unc_diff']:12.4f} {val_str:>10s} {ratio_str:>10s}")
    return rows


    # plot uncertainty differences
    fig, ax = plt.subplots(1, 1, figsize=(max(6, len(names) * 0.6), max(4, len(all_vals) * 0.2)))
    im = ax.imshow(Zu, aspect='auto', origin='lower', interpolation='none')
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90)
    ax.set_yticks(range(len(all_vals)))
    ax.set_yticklabels([f'{v:.3f}' for v in all_vals])
    ax.set_xlabel('Scan')
    ax.set_ylabel('Correlation')
    ax.set_title('Total unc. - nominal unc. [GeV]')
    fig.colorbar(im, ax=ax)
    # avoid tight_layout warning for mixed Axes; use subplots_adjust and expand width
    fig.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.12)
    p2_png = os.path.join(outdir, '2d_uncertainty_diff.png')
    p2_pdf = os.path.join(outdir, '2d_uncertainty_diff.pdf')
    fig.savefig(p2_png)
    try:
        fig.savefig(p2_pdf)
    except Exception:
        pass
    plt.close(fig)



# a function that takes the json and then creates a latex table sumamrizing for each correlation (one per row) the max central positive and negative deviation and max pos. and neg. uncertainty deviation from the central value
def create_latex_summary_table(scans: Dict[str, dict], outdir: str, unblind: bool = False) -> None:
    """Create a LaTeX summary table of max positive/negative central and uncertainty deviations per scan.
    Sorted by impact on central value, then by impact on uncertainty.
    The table is written to 'scan_summary_table.tex' in outdir.
    """
    table_lines = []


    table_lines.append(r'\begin{table}[h!]')
    table_lines.append(r'\scriptsize')
    table_lines.append(r'\caption{Summary of maximum positive and negative deviations in fitted mass and uncertainty for each correlation scan.}')
    table_lines.append(r'\resizebox{\textwidth}{!}{%')

    table_lines.append(r'\begin{tabular}{lcccc}')
    table_lines.append(r'\hline')
    table_lines.append(r'Scan & Max $+\Delta m_{t}$ [GeV] & Max $-\Delta m_{t}$ [GeV] & Max $+\Delta \sigma$ [GeV] & Max $-\Delta \sigma$ [GeV] \\')
    table_lines.append(r'\hline')




    # Sort scans by max impact on central value, then by max impact on uncertainty
    def sort_key(name):
        s = scans[name]
        if not isinstance(s, dict):
            return (0.0, 0.0)
        entries = s.get('entries', []) or []

        max_pos_delta = float('-inf')
        max_neg_delta = float('inf')
        max_pos_unc_diff = float('-inf')
        max_neg_unc_diff = float('inf')

        for e in entries:
            c = e.get('central')
            u = e.get('total_unc')
            nominal = s.get('nominal')
            nom_c = nominal.get('central') if nominal else None
            nom_u = nominal.get('total_unc') if nominal else None

            if c is not None and nom_c is not None:
                delta = c - nom_c
                if delta > max_pos_delta:
                    max_pos_delta = delta
                if delta < max_neg_delta:
                    max_neg_delta = delta

            if u is not None and nom_u is not None:
                unc_diff = u - nom_u
                if unc_diff > max_pos_unc_diff:
                    max_pos_unc_diff = unc_diff
                if unc_diff < max_neg_unc_diff:
                    max_neg_unc_diff = unc_diff

        return (max(abs(max_pos_delta), abs(max_neg_delta)), max(abs(max_pos_unc_diff), abs(max_neg_unc_diff)))
    
    sorted_names = sorted(scans.keys(), key=sort_key, reverse=True)




    for name in sorted_names:
        if "fit.log" in name.lower() or name.lower() == "nominal":
            continue
        s = scans[name]
        if not isinstance(s, dict):
            continue
        entries = s.get('entries', []) or []

        max_pos_delta = float('-inf')
        max_neg_delta = float('inf')
        max_pos_unc_diff = float('-inf')
        max_neg_unc_diff = float('inf')

        for e in entries:
            c = e.get('central')
            u = e.get('total_unc')
            nominal = s.get('nominal')
            nom_c = nominal.get('central') if nominal else None
            nom_u = nominal.get('total_unc') if nominal else None

            if c is not None and nom_c is not None:
                delta = c - nom_c
                if delta > max_pos_delta:
                    max_pos_delta = delta
                if delta < max_neg_delta:
                    max_neg_delta = delta

            if u is not None and nom_u is not None:
                unc_diff = u - nom_u
                if unc_diff > max_pos_unc_diff:
                    max_pos_unc_diff = unc_diff
                if unc_diff < max_neg_unc_diff:
                    max_neg_unc_diff = unc_diff

        # format values for LaTeX
        pos_delta_str = f"{max_pos_delta:.3f}" if max_pos_delta != float('-inf') else 'n/a'
        neg_delta_str = f"{max_neg_delta:.3f}" if max_neg_delta != float('inf') else 'n/a'
        pos_unc_str = f"{max_pos_unc_diff:.3f}" if max_pos_unc_diff != float('-inf') else 'n/a'
        neg_unc_str = f"{max_neg_unc_diff:.3f}" if max_neg_unc_diff != float('inf') else 'n/a'
        # take care of backspacing underscores in name
        name_escaped = name.replace('_', r'\_')
        # now let's use three different colors: black CMS, red ATLAS 8 TeV and blue ATLAS 13 TeV
        # we split each name in the middle with "__" and then color both sides accordingly
        name_colored = ''
        if '\_\_' in name_escaped:
            parts = name_escaped.split('\_\_', 1)
            if len(parts) == 2:
                left, right = parts
                # print ('DEBUG: left =', left, ' right =', right)
                if 'ATLAS13' in right or 'ATLAS\_13' in right or 'ATLAS13TeV' in right:
                    right = f"{{\\color{{blue}} {right}}}"
                if 'ATLAS13' in left or 'ATLAS\_13' in left or 'ATLAS13TeV' in left:
                    left = f"{{\\color{{blue}} {right}}}"
                if 'ATLAS8' in right or 'ATLAS\_8' in right or 'ATLAS8TeV' in right:
                    right = f"{{\\color{{red}} {right}}}"
                if 'ATLAS8' in left or 'ATLAS\_8' in left or 'ATLAS8TeV' in left:
                    left = f"{{\\color{{red}} {left}}}"
                name_colored = f"{left} {right}"
            else:
                name_colored = name_escaped
        else:
            name_colored = name_escaped

        table_lines.append(f"{name_colored} & {pos_delta_str} & {neg_delta_str} & {pos_unc_str} & {neg_unc_str} \\\\")
    table_lines.append(r'\hline')
    table_lines.append(r'\end{tabular}')
    table_lines.append(r'}')  # end resizebox
    table_lines.append(r'\end{table}')
    os.makedirs(outdir, exist_ok=True)
    table_path = os.path.join(outdir, 'scan_summary_table.tex')
    with open(table_path, 'w') as tf:
        tf.write('\n'.join(table_lines))
    if os.path.isfile(table_path):
        print('Wrote LaTeX summary table to', table_path)
    


def main():
    parser = argparse.ArgumentParser(description='Collect Condor scan fit.log results and plot summaries')
    parser.add_argument('jobs_folder', help='Jobs folder created by createCondorJobs.py')
    parser.add_argument('--eos-output', help='Optional EOS output base to also scan (where jobs copied outputs)')
    parser.add_argument('--outdir', default='condor_scan_plots', help='Directory where plots and JSON are written')
    parser.add_argument('--json', default='scan_summary.json', help='JSON output filename (inside outdir)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose messages')
    parser.add_argument('-d', '--debug', action='store_true', help='Enable debug output')
    parser.add_argument('-1d', '--plot1d', action='store_true', help='Enable 1D plots')
    parser.add_argument('-2d', '--plot2d', action='store_true', help='Enable 2D plots')
    parser.add_argument('-tex', '--maketex', action='store_true', help='Enable LaTeX table output')
    parser.add_argument('--load-json', help='Path to existing scan JSON to load (skip parsing if exists)')
    parser.add_argument('--unblind', action='store_true', default = False, help='Unblind plots (show actual mass values instead of deltas)')
    args = parser.parse_args()

    scans = None
    if args.load_json:
        load_path = os.path.abspath(args.load_json)
        if os.path.isfile(load_path):
            if args.verbose:
                print(f'Loading scans from JSON: {load_path}')
            with open(load_path, 'r') as lj:
                scans = json.load(lj)
        else:
            if args.verbose:
                print(f'Warning: --load-json provided but file not found: {load_path} -- will parse results')

    if scans is None:
        scans = collect_results(args.jobs_folder, args.eos_output, verbose=args.verbose, debug=args.debug)

    # normalize structure: promote fit.log/value=None entries to nominal and ensure nominal in entries
    try:
        normalize_scans(scans, verbose=args.verbose)
    except Exception:
        if args.debug:
            print('DEBUG: normalize_scans failed')
    # propagate a global nominal scan into individual scans if available
    try:
        propagate_global_nominal(scans, verbose=args.verbose)
    except Exception:
        if args.debug:
            print('DEBUG: propagate_global_nominal failed')

    # dump names mapping for 2D plots (short labels)
    try:
        print ("Writing names mapping for 2D plots...")
        names = sorted(scans.keys())
        namesNew = [str(i) for i in range(len(names))]
        mapping_path = os.path.join(args.outdir, 'names_mapping.txt')
        with open(mapping_path, 'w') as mf:
            for i, n in enumerate(names):
                mf.write(f"{i}\t{n}\n")
        if args.verbose:
            print('Wrote names mapping to', mapping_path)
    except Exception:
        if args.debug:
            import traceback
            traceback.print_exc()

    # compute per-scan linear correlations and write ranking
    try:
        corr_rows = compute_scan_correlations(scans, verbose=args.verbose)
        # print top 10
        print('\nTop 10 scans by absolute Pearson r (mass vs correlation):')
        print('{:3s} {:60s} {:>8s} {:>10s} {:>10s} {:>10s}'.format('#', 'scan', 'r', 'slope', 'intercept', 'n'))
        for i, r in enumerate(corr_rows[:10]):
            print(f"{i+1:3d} {r['scan'][:60]:60s} {r['r']:8.4f} {r['slope']:10.4f} {r['intercept']:10.4f} {r['n_points']:10d}")
        # write CSV
        import csv
        cfn = os.path.join(args.outdir, 'rank_mass_correlation.csv')
        with open(cfn, 'w', newline='') as cf:
            writer = csv.DictWriter(cf, fieldnames=['scan', 'r', 'slope', 'intercept', 'n_points'])
            writer.writeheader()
            for r in corr_rows:
                writer.writerow(r)
        if args.verbose:
            print('Wrote mass-correlation ranking to', cfn)
    except Exception:
        if args.debug:
            import traceback
            traceback.print_exc()

    os.makedirs(args.outdir, exist_ok=True)
    jfn = os.path.join(args.outdir, args.json)
    with open(jfn, 'w') as jf:
        json.dump(scans, jf, indent=2)


    create_latex_summary_table(scans, args.outdir, unblind=args.unblind)


    # per-scan plots
    if args.plot1d:
        # clean outputs first
        for f in os.listdir(args.outdir):
            if f.startswith('scan_') and (f.endswith('.png') or f.endswith('.pdf')):
                try:
                    os.remove(os.path.join(args.outdir, f))
                except Exception:
                    pass
        for name, scan in scans.items():
            plot_scan(name, scan, args.outdir, unblind=args.unblind)

    # 2D summaries
    if args.plot2d:
        # clean
        for f in os.listdir(args.outdir):
            if f.startswith('2d_') and (f.endswith('.png') or f.endswith('.pdf')):
                try:
                    os.remove(os.path.join(args.outdir, f))
                except Exception:
                    pass
        plot_2d(scans, args.outdir)

    # print top deviations summary
    try:
        central_rows = print_top_deviations(scans, topn=10, unblind=args.unblind)
    except Exception:
        if args.debug:
            print('DEBUG: print_top_deviations failed')

    try:
        uncertainty_rows = print_top_uncertainty_deviations(scans, topn=10)
    except Exception:
        if args.debug:
            print('DEBUG: print_top_uncertainty_deviations failed')

    # write full rankings to CSV
    try:
        import csv
        cfn = os.path.join(args.outdir, 'rank_central_deviation.csv')
        with open(cfn, 'w', newline='') as cf:
            writer = csv.DictWriter(cf, fieldnames=['scan', 'max_abs_delta', 'value', 'path', 'ratio_nominal_unc'])
            writer.writeheader()
            for r in (central_rows or []):
                writer.writerow({
                    'scan': r.get('scan'),
                    'max_abs_delta': r.get('max_abs_delta'),
                    'value': r.get('value'),
                    'path': r.get('path'),
                    'ratio_nominal_unc': r.get('ratio_nominal_unc'),
                })

        ufn = os.path.join(args.outdir, 'rank_uncertainty_deviation.csv')
        with open(ufn, 'w', newline='') as uf:
            writer = csv.DictWriter(uf, fieldnames=['scan', 'max_abs_unc_diff', 'value', 'path', 'ratio_nominal_unc'])
            writer.writeheader()
            for r in (uncertainty_rows or []):
                writer.writerow({
                    'scan': r.get('scan'),
                    'max_abs_unc_diff': r.get('max_abs_unc_diff'),
                    'value': r.get('value'),
                    'path': r.get('path'),
                    'ratio_nominal_unc': r.get('ratio_nominal_unc'),
                })
        if args.verbose:
            print('Wrote ranking CSVs to', os.path.abspath(args.outdir))
    except Exception:
        if args.debug:
            import traceback
            traceback.print_exc()

    print('Wrote JSON and plots to', os.path.abspath(args.outdir))


if __name__ == '__main__':
    main()
