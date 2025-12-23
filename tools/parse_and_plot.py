#!/usr/bin/env python3
"""
Parser and plotting tool for Convino-style input files.

Usage:
    python tools/parse_and_plot.py <input_file> [--outdir plots]

Produces one PNG per uncertainty (column) in <outdir>.

Behavior:
- Parses sections like [estimates] and [not fitted].
- From [estimates] reads name_i / value_i pairs and uses x positions 0..N-1.
- From [not fitted] expects a header line listing column (uncertainty) names, then rows where
  each row starts with a central-name followed by numeric shifts for each column.
- For each uncertainty column, a figure is created with two panels:
    top: central values (black) and shifted values (red)
    bottom: ratio = shift / central (red dots)

This script is intentionally simple and tolerant to whitespace.
"""
import argparse
import os
import sys
import math
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep

# hep.style.use('CMS')

def parse_sections(lines):
    sections = OrderedDict()
    current = None
    buf = []
    for raw in lines:
        line = raw.rstrip('\n')
        if line.strip().startswith('[') and line.strip().endswith(']'):
            if current is not None:
                sections[current] = buf
            current = line.strip()[1:-1].strip()
            buf = []
        else:
            if current is not None:
                buf.append(line)
    if current is not None:
        sections[current] = buf
    return sections


def parse_estimates(lines):
    # lines contains the lines between [estimates] and [end estimates]
    names = []
    values = []
    for l in lines:
        s = l.strip()
        if not s:
            continue
        if s.startswith('name_'):
            # format: name_i =  ATLAS_... (may have spacing)
            parts = s.split('=', 1)
            if len(parts) == 2:
                name = parts[1].strip()
                names.append(name)
        elif s.startswith('value_'):
            parts = s.split('=', 1)
            if len(parts) == 2:
                try:
                    val = float(parts[1].strip())
                except Exception:
                    val = float('nan')
                values.append(val)
    # If values length equals names length, build ordered mapping
    if len(names) != len(values):
        # fallback: try to parse paired blocks sequentially
        paired = []
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith('name_'):
                name = line.split('=',1)[1].strip()
                # look ahead for value
                j = i+1
                val = float('nan')
                while j < len(lines):
                    if lines[j].strip().startswith('value_'):
                        val = float(lines[j].split('=',1)[1].strip())
                        break
                    j += 1
                paired.append((name, val))
                i = j+1
            else:
                i += 1
        if paired:
            names, values = zip(*paired)
            return list(names), list(values)
    return names, values


def parse_notfitted(lines):
    # Expect first non-empty line to be header of column names (systematics)
    header_tokens = []
    data_rows = []  # tuples (row_name, [numbers...])
    i = 0
    # skip initial blank lines
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines):
        return header_tokens, data_rows
    # header may span multiple lines if long; but we will assume it's the first non-empty line
    header_line = lines[i].strip()
    header_tokens = header_line.split()
    i += 1
    # remaining lines until end: rows, each starting with a row name and then numeric columns
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line:
            continue
        if line.startswith('#'):
            continue
        parts = line.split()
        row_name = parts[0]
        nums = []
        for p in parts[1:]:
            try:
                nums.append(float(p))
            except Exception:
                # ignore non-numeric tokens
                pass
        data_rows.append((row_name, nums))
    return header_tokens, data_rows


def build_uncertainty_dict(names_order, header_tokens, data_rows):
    # names_order: list of central names in order (from estimates)
    # header_tokens: list of column names
    # data_rows: list of (row_name, [nums...]) where row_name corresponds to a central name

    # Build mapping row_name->nums
    rowmap = {rname: nums for (rname, nums) in data_rows}

    # Number of columns
    if len(header_tokens) == 0:
        # nothing to do
        return OrderedDict()
    ncols = len(header_tokens)

    # Prepare uncertainties dict: key -> list of values aligned with names_order
    unc = OrderedDict()
    for col in header_tokens:
        unc[col] = []

    for nm in names_order:
        nums = rowmap.get(nm)
        if nums is None:
            # If exact name not found try to match with prefix
            found = None
            for k in rowmap.keys():
                if k.startswith(nm):
                    found = k
                    break
            if found:
                nums = rowmap[found]
            else:
                # fill with NaNs
                nums = [float('nan')] * ncols
        # If row has fewer columns than header, pad
        if len(nums) < ncols:
            nums = nums + [float('nan')] * (ncols - len(nums))
        for j, col in enumerate(header_tokens):
            unc[col].append(nums[j])
    return unc


def plot_uncertainty(x, names, central_values, shifts, unc_name, outpath):
    # x: list or array positions
    # names: list labels
    # central_values: list of central floats
    # shifts: list of shift floats (same length)

    central = np.array(central_values, dtype=float)
    shifts = np.array(shifts, dtype=float)
    shifted_values = central + shifts

    fig = plt.figure(figsize=(10,6))
    gs = fig.add_gridspec(2,1, height_ratios=[3,1], hspace=0.05)
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1], sharex=ax_top)

    # top: central (black points) and shifted (red points)
    ax_top.plot(x, central, 'k.-', label='central')
    ax_top.plot(x, shifted_values, 'ro', label=f'shifted ({unc_name})')
    # connect with vertical lines
    for xi, c, s in zip(x, central, shifted_values):
        ax_top.vlines(xi, min(c,s), max(c,s), colors='red', alpha=0.4)

    ax_top.set_ylabel('Value')
    ax_top.legend()
    ax_top.grid(True)

    # bottom: ratio = shift / central (handle zeros)
    ratio = np.full_like(central, np.nan, dtype=float)
    for i in range(len(central)):
        c = central[i]
        if c != 0 and not math.isnan(c):
            ratio[i] = (shifted_values[i]) / c
        else:
            ratio[i] = float('nan')

    ax_bot.axhline(0, color='gray', linewidth=0.8)
    ax_bot.plot(x, ratio, 'ro-')
    ax_bot.set_ylabel('shift/central')
    ax_bot.set_xlabel('index')
    ax_bot.grid(True)

    # set xticks to names (rotate)
    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(names, rotation=45, ha='right', fontsize=8)

    # center around 1.0 and set range according to uncertainty size
    valid_ratios = ratio[~np.isnan(ratio)]
    if len(valid_ratios) > 0:
        rmin = np.min(valid_ratios)
        rmax = np.max(valid_ratios)
        rcenter = 1.0
        rspan = max(abs(rmax - rcenter), abs(rcenter - rmin))
        margin = rspan * 0.2
        ax_bot.set_ylim(rcenter - rspan - margin, rcenter + rspan + margin)
    else:
        ax_bot.set_ylim(0.8, 1.2)

    fig.suptitle(f'Uncertainty: {unc_name}')
    # plt.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(outpath)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Parse Convino input and plot uncertainties per-column')
    parser.add_argument('input_file', help='Path to input file')
    parser.add_argument('--outdir', default='plots', help='Output directory for plots')
    parser.add_argument('--max-plots', type=int, default=None, help='Limit number of plots (useful for large files)')
    args = parser.parse_args()

    # check if outdir exists, if so, clean it, and delete plots inside
    if os.path.exists(args.outdir):
        for f in os.listdir(args.outdir):
            if f.endswith('.png'):
                os.remove(os.path.join(args.outdir, f))
    else:
        os.makedirs(args.outdir)

    with open(args.input_file, 'r') as f:
        lines = f.readlines()
    sections = parse_sections(lines)

    # Parse estimates
    est_lines = sections.get('estimates', [])
    names, values = parse_estimates(est_lines)
    if not names:
        print('No estimates parsed; aborting')
        sys.exit(1)

    # Parse not fitted
    notf_lines = sections.get('not fitted', [])
    header_tokens, data_rows = parse_notfitted(notf_lines)

    unc = build_uncertainty_dict(names, header_tokens, data_rows)

    # from the dict make a ranking of uncertainties by max absolute shift
    unc_ranking = sorted(unc.items(), key=lambda item: max(abs(s) for s in item[1] if not math.isnan(s)), reverse=True)
    unc = OrderedDict(unc_ranking)
    # print summary, sorted by max shift, but printas relative to central values
    print('Uncertainty summary (sorted by max absolute shift):')
    for unc_name, shifts in unc.items():
        max_rel = 0.0
        for i in range(len(shifts)):
            c = values[i]
            s = shifts[i]
            if c != 0 and not math.isnan(c) and not math.isnan(s):
                rel = abs(s) / abs(c)
                if rel > max_rel:
                    max_rel = rel
        print(f'  {unc_name}: max relative shift = {max_rel*100:.2f} %')



    # Prepare output
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    x = list(range(len(names)))

    # Iterate uncertainties
    count = 0
    for unc_name, shifts in unc.items():
        if args.max_plots is not None and count >= args.max_plots:
            break
        # outpath = os.path.join(outdir, f"unc_{count:03d}_{unc_name}.png")
        outpath = os.path.join(outdir, f"unc_{unc_name}.png")
        print(f'Plotting {unc_name} -> {outpath} (n={len(shifts)})')
        plot_uncertainty(x, names, values, shifts, unc_name, outpath)
        count += 1

    print(f'Done. Plots written to {outdir}')

if __name__ == '__main__':
    main()
