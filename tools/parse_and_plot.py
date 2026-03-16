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

hep.style.use('CMS')

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
    # we need to leave more space below the bottom axes plot for the rotated x labels
    plt.subplots_adjust(bottom=0.2)
    gs = fig.add_gridspec(2,1, height_ratios=[3,1], hspace=0.05)
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[1], sharex=ax_top)

    # top: central (black points) and shifted (red points)
    ax_top.plot(x, central, 'k.-', label='Central')
    ax_top.plot(x, shifted_values, 'ro', label=f'{unc_name}')
    # connect with vertical lines
    for xi, c, s in zip(x, central, shifted_values):
        ax_top.vlines(xi, min(c,s), max(c,s), colors='red', alpha=0.4)

    # no x labels for the top axes
    plt.setp(ax_top.get_xticklabels(), visible=False)

    ax_top.set_ylabel('Cross section')
    # set y range large enough
    all_vals = np.concatenate([central, shifted_values])
    vmin = np.nanmin(all_vals)
    vmax = np.nanmax(all_vals)
    vrange = vmax - vmin
    ax_top.set_ylim(0., vmax + 0.1*vrange)
    ax_top.legend(fontsize=9, loc ='upper right')
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
    # ax_bot.set_xlabel('index')
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

    hep.cms.label(exp="ATLAS+CMS", llabel= "Work in Progress", rlabel = "8+13 TeV", ax=ax_top)
    # fig.suptitle(f'Uncertainty: {unc_name}')
    # plt.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(outpath)
    # also save as pdf
    pdf_outpath = outpath.rsplit('.',1)[0] + '.pdf'
    fig.savefig(pdf_outpath)
    plt.close(fig)


# define a function that plots all uncertainties in one single plot, sorted by max absolute shift summed
def plot_all_uncertainties_relative(x, names, central_values, unc_dict, outpath):
    fig, ax = plt.subplots(figsize=(12,8))

    central = np.array(central_values, dtype=float)

    # ensure that I have enough colors and plottung styles for the length of the group
    colors = plt.cm.get_cmap('tab10').colors
    linestyles = ['-', '--', '-.', ':']
    color_cycle = iter(colors)
    linestyle_cycle = iter(linestyles)

    for unc_name, shifts in unc_dict.items():
        shifts = np.array(shifts, dtype=float)
        shifted_values = central + shifts
        # compute max absolute shift
        max_abs_shift = np.max(np.abs(shifts/central[~np.isnan(shifts)]))
        ax.plot(x, shifted_values/central, label=f'{unc_name} (max rel. shift: {max_abs_shift:.2f})', alpha=1.0, color=next(color_cycle))

    # ax.plot(x, central, 'k.-', label='central', linewidth=2)

    ax.set_ylabel('Relative uncertainties')
    # ax.set_xlabel('Bins')
    # ax.set_title('All Uncertainties')
    # ax.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.15, 1))
    ax.legend(fontsize=9, loc='upper left', ncol = 2)
    ax.grid(True)
    # center around 1.0
    ax.axhline(1.0, color='gray', linewidth=0.8)
    max_shift_vals = max([np.max(np.abs((shifts+central)/central)) for shifts in unc_dict.values()])
    min_shift_vals = min([np.min(np.abs((shifts+central)/central)) for shifts in unc_dict.values()])
    ax.set_ylim(0.95*min_shift_vals, 1.05*max_shift_vals)

    # set xticks to names (rotate)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=9)

    hep.cms.label(exp="ATLAS+CMS", llabel= "Work in Progress", rlabel = "8+13 TeV", ax=ax)

    plt.tight_layout()
    fig.savefig(outpath)
    # also save as pdf
    pdf_outpath = outpath.rsplit('.',1)[0] + '.pdf'
    fig.savefig(pdf_outpath)
    plt.close(fig)




def main():
    parser = argparse.ArgumentParser(description='Parse Convino input and plot uncertainties per-column')
    parser.add_argument('input_file', help='Path to input file')
    parser.add_argument('--outdir', default='plots', help='Output directory for plots')
    parser.add_argument('--max-plots', type=int, default=None, help='Limit number of plots (useful for large files)')
    parser.add_argument('--plot-all', action='store_true', help='Plot all uncertainties in one plot')
    parser.add_argument('--plot-single', action='store_true', help='Plot single uncertainties in individual plots')
    args = parser.parse_args()

    if not args.plot_all and not args.plot_single:
        print('No plotting option specified (--plot-all or --plot-single); aborting')
        sys.exit(1)

    # check if outdir exists, if so, clean it, and delete plots inside
    if os.path.exists(args.outdir):
        for f in os.listdir(args.outdir):
            if f.endswith('.png') or f.endswith('.pdf') or f.endswith('.csv'):
                os.remove(os.path.join(args.outdir, f))
        if args.plot_all:
            grouped_dir = os.path.join(args.outdir, 'grouped')
            if os.path.exists(grouped_dir):
                for f in os.listdir(grouped_dir):
                    if f.endswith('.png') or f.endswith('.pdf'):
                        os.remove(os.path.join(grouped_dir, f))
            else:
                os.makedirs(grouped_dir)
        if args.plot_single:
            individual_dir = os.path.join(args.outdir, 'individual')
            if os.path.exists(individual_dir):
                for f in os.listdir(individual_dir):
                    if f.endswith('.png') or f.endswith('.pdf'):
                        os.remove(os.path.join(individual_dir, f))
            else:
                os.makedirs(individual_dir)
    else:
        os.makedirs(args.outdir)
        if args.plot_all:
            os.makedirs(os.path.join(args.outdir, 'grouped'), exist_ok=True)
        if args.plot_single:
            os.makedirs(os.path.join(args.outdir, 'individual'), exist_ok=True)

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


    # write an output csv file with the absolute and relative shifts for ALL bins in one entry each, sorted by overall summed magnitude. One bin of the distributuon is one column. In addition we have the sum columns.
    csv_outpath = os.path.join(args.outdir, 'uncertainty_summary.csv')
    with open(csv_outpath, 'w') as cf:
        # header
        cf.write('uncertainty_name')
        for i in range(len(names)):
            cf.write(f',bin_{i}_abs,bin_{i}_rel')
        cf.write(',sum_abs,sum_rel\n')
        # rows
        for unc_name, shifts in unc.items():
            cf.write(unc_name)
            sum_abs = 0.0
            sum_rel = 0.0
            for i in range(len(shifts)):
                s = shifts[i]
                c = values[i]
                abs_s = s if not math.isnan(s) else 0.0
                rel_s = (abs(s) / abs(c)) if (c != 0 and not math.isnan(c) and not math.isnan(s)) else 0.0
                sum_abs += abs(abs_s)
                sum_rel += rel_s
                cf.write(f',{abs_s},{rel_s}')
            cf.write(f',{sum_abs},{sum_rel}\n')
    print(f'Wrote uncertainty summary to {csv_outpath}')


    # Prepare output
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    x = list(range(len(names)))

    # Iterate uncertainties
    if args.plot_all:
        outpath = os.path.join(outdir, 'all_uncertainties.png')
        print(f'Plotting all uncertainties -> {outpath}')
        # split set of uncertainties in groups of 10 each
        unc_groups = [list(unc.items())[i:i + 10] for i in range(0, len(unc), 10)]
        for igroup, group in enumerate(unc_groups):
            group_dict = OrderedDict(group)
            outpath_group = os.path.join(outdir, f'grouped/all_uncertainties_group_{igroup:02d}.png')
            print(f' Plotting group {igroup+1}/{len(unc_groups)} -> {outpath_group}')
            plot_all_uncertainties_relative(x, names, values, group_dict, outpath_group)
    if args.plot_single:
        count = 0
        for unc_name, shifts in unc.items():
            if args.max_plots is not None and count >= args.max_plots:
                break
            # outpath = os.path.join(outdir, f"unc_{count:03d}_{unc_name}.png")
            outpath = os.path.join(outdir, f"individual/unc_{unc_name}.png")
            print(f'Plotting {unc_name} -> {outpath} (n={len(shifts)})')
            plot_uncertainty(x, names, values, shifts, unc_name, outpath)
            count += 1

    print(f'Done. Plots written to {outdir}')

if __name__ == '__main__':
    main()
