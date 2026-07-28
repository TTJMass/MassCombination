#!/usr/bin/env python3
"""Generate the v1 order/PDF/dataset/POI-compatibility fit matrix as Condor jobs.

Sits above createCondorJobs.py (left unmodified -- it's production
infrastructure for the correlation-scan axis) rather than extending its core
loop: for each axis-tuple in the concrete v1 matrix (filtered through
matrix_axes.valid_combo), calls createCondorJobs.py --only-nominal once with
--dofit-extra-args assembled to select dataset/PDF/order/POI-config.

Two things createCondorJobs.py does NOT vary on its own that matter here:
  - --eos-output-path: createCondorJobs.py derives its EOS output dir purely
    from the setup directory's basename, with no awareness of the PDF/order/
    poi_config axis. Calling it multiple times over the *same* setup (e.g.
    the 15-job PDF/order sweep) with a shared eos-output-path would silently
    collide -- every job's fit.log/.completed marker would land in the same
    EOS folder and overwrite the previous one. This script gives every
    axis-tuple its own eos_output_path, not just its own jobs_folder.
  - --exp: the real dataset_key (not createCondorJobs.py's synthetic
    batch_name/jobname label) is passed via --dofit-extra-args, appended
    after the SH template's own hardcoded --exp so argparse's last-flag-wins
    behavior lets it take over. This requires fit_object.py's expPath-before-
    input_f-lookup precedence fix (see fit_object.py __init__) -- without it,
    a real dataset_key --exp would make the fit silently ignore --expPath
    (the job's freshly-produced convino result) and try to read a stale,
    nonexistent path on the worker instead.
"""
import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'mtpole-ttj-pyconvino'))
from matrix_axes import (SETUP_TO_DATASET_KEY, DATASET_KEY_TO_SETUP, theory_tag, valid_combo,  # noqa: E402
                          REF_DATASET_KEY, REF_THEORY_SOURCE, REF_ORDER, REF_PDF,
                          SWEEP_DATASET_KEYS, SWEEP_ORDERS, SWEEP_PDFS, POLY_ORDERS)
from createCondorJobs import sanitize  # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

STANDALONE_DATASET_KEYS = ['CMS_13TeV_npz', 'ATLAS_8TeV_npz', 'ATLAS_13TeV_npz']
# ATLAS-only energy combination (8+13 TeV, no CMS) -- previously had NO
# single (non-split) legacy-NLO entry anywhere in the matrix at all; only
# existed via POI_SPLIT_ENTRIES's 'split_meas' fit, a different fit. Added
# 2026-07-13 so it has a legacy point like every other one of the 7 user-facing
# fits (needed for its own diagnostic plots and its own legacy-vs-new-NLO row).
ATLAS_ENERGY_COMBO_DATASET_KEY = 'ATLAS_813TeV_npz'
OTHER_COMBO_DATASET_KEYS = ['ATLAS_13TeV_CMS_13TeV_npz', 'ATLAS_8TeV_CMS_13TeV_npz']

STRIPPER_ORDERS_FOR_SWEEP = SWEEP_ORDERS
STRIPPER_PDFS_FOR_SWEEP = SWEEP_PDFS

# (dataset_key, poi_config flag) for the v1 POI-compatibility matrix -- all
# at the reference NLO/CT18NLO point only (plan decision #6).
POI_SPLIT_ENTRIES = [
    # CMS-vs-ATLAS (2-POI): only --splitMassesExp separates ATLAS from CMS.
    (REF_DATASET_KEY, 'split_exp'),
    # ATLAS8-vs-ATLAS13 (2-POI): single experiment, two energies -- only
    # --splitMasses (per-measurement) separates them; --splitMassesExp would
    # collapse both into one 'ATLAS' mass parameter.
    ('ATLAS_813TeV_npz', 'split_meas'),
    # ATLAS8-vs-CMS13 (2-POI).
    ('ATLAS_8TeV_CMS_13TeV_npz', 'split_meas'),
    # Full 3-POI split.
    (REF_DATASET_KEY, 'split_meas'),
]


def build_matrix():
    """Returns a list of axis-tuple dicts making up the concrete v1 matrix."""
    entries = []

    def add(category, dataset_key, theory_source, order, pdf, poi_config):
        entries.append({
            'category': category, 'dataset_key': dataset_key,
            'theory_source': theory_source, 'order': order, 'pdf': pdf,
            'poi_config': poi_config,
        })

    add('reference', REF_DATASET_KEY, REF_THEORY_SOURCE, REF_ORDER, REF_PDF, 'single')

    for dk in STANDALONE_DATASET_KEYS:
        add('standalone', dk, REF_THEORY_SOURCE, REF_ORDER, REF_PDF, 'single')

    add('atlas_energy_combo', ATLAS_ENERGY_COMBO_DATASET_KEY, REF_THEORY_SOURCE, REF_ORDER, REF_PDF, 'single')

    for dk in OTHER_COMBO_DATASET_KEYS:
        add('other_combo', dk, REF_THEORY_SOURCE, REF_ORDER, REF_PDF, 'single')

    # Full stripper order x PDF sweep for all 7 user-facing fits (extended
    # 2026-07-13 from just the reference/full-combo dataset) -- drives each
    # fit's own order-to-order / PDF-to-PDF / legacy-vs-new-NLO diagnostic
    # plots in matrix_comparison_plot.py.
    for dk in SWEEP_DATASET_KEYS:
        for order in STRIPPER_ORDERS_FOR_SWEEP:
            for pdf in STRIPPER_PDFS_FOR_SWEEP:
                add('sweep', dk, 'stripper_json', order, pdf, 'single')

    for dk, poi_config in POI_SPLIT_ENTRIES:
        add('poi_split', dk, REF_THEORY_SOURCE, REF_ORDER, REF_PDF, poi_config)

    return entries


def axis_tag(entry):
    setup = DATASET_KEY_TO_SETUP[entry['dataset_key']]
    return '__'.join([setup, theory_tag(entry['theory_source'], entry['order']), entry['pdf'], entry['poi_config'],
                       f"poly{entry['poly_order']}"])


def dofit_extra_args_for(entry, common_extra_args):
    parts = [common_extra_args] if common_extra_args else []
    parts.append(f"--exp {entry['dataset_key']}")
    parts.append(f"--PDF {entry['pdf']}")
    if entry['theory_source'] == 'stripper_json':
        parts.append('--stripper')
        parts.append('--massDepVariations')
    parts.append(f"--order {entry['order']}")
    if entry['poi_config'] == 'split_meas':
        parts.append('--splitMasses')
    elif entry['poi_config'] == 'split_exp':
        parts.append('--splitMassesExp')
    parts.append(f"--polyOrder {entry['poly_order']}")
    return ' '.join(parts)


def main():
    parser = argparse.ArgumentParser(description='Generate the v1 fit matrix as Condor jobs (drives createCondorJobs.py per axis-tuple)')
    parser.add_argument('eos_output_base', help='Base EOS path; each axis-tuple gets its own subdirectory under here')
    parser.add_argument('jobs_folder_base', help='Base local folder; each axis-tuple gets its own jobs_folder subdirectory under here')
    parser.add_argument('--batch-name', default=None,
                         help='Condor batch_name prefix (default: MassComb_fitmatrix_<YYYYMMDD>)')
    parser.add_argument('--categories', nargs='+',
                         choices=['reference', 'standalone', 'atlas_energy_combo', 'sweep', 'other_combo',
                                  'poi_split', 'all'],
                         default=['all'], help='Restrict generation to these categories')
    parser.add_argument('--common-dofit-extra-args', default='--interpCheck --nuisanceFit --pulls --budget',
                         help='Extra doFit.py args applied to every job, before the per-axis-tuple ones '
                              '(mirrors writeAllCondorJobs_pyconvino.sh commonArgs)')
    parser.add_argument('--conda-pack-tarball',
                         default=os.path.join(REPO_ROOT, 'envCache', 'masscomb_packed.tar.gz'),
                         help='Passed through to createCondorJobs.py')
    parser.add_argument('--nlo-theory-json',
                         default=os.path.join(REPO_ROOT, 'theory-data', 'nlo_converted.json'),
                         help='Passed through to createCondorJobs.py')
    parser.add_argument('--stripper-theory-json',
                         default=os.path.join(REPO_ROOT, 'theory-data', 'newdata.json'),
                         help='Passed through to createCondorJobs.py')
    parser.add_argument('--blind-salt-file', default=os.path.expanduser('~/.masscomb_blind_salt'),
                         help='Passed through to createCondorJobs.py')
    parser.add_argument('--do-impacts', dest='do_impacts', action='store_true', default=True,
                         help='Passed through to createCondorJobs.py (default: on, mirrors '
                              'writeAllCondorJobs_pyconvino.sh commonArgs -- the fit matrix is O(25) '
                              'jobs, not the O(5000)-job correlation scans --no-impacts exists for, '
                              'and every job needs a real stat_only_covariance: --no-impacts makes '
                              'pyconvino write NaN placeholders there, which crashes exp_xsec.py '
                              'downstream -- see condor_fitmatrix_eigh_crash memory)')
    parser.add_argument('--no-do-impacts', dest='do_impacts', action='store_false',
                         help='Escape hatch: fall back to createCondorJobs.py\'s --no-impacts default')
    parser.add_argument('--dry-run', action='store_true', help='Print the resolved matrix and commands without invoking createCondorJobs.py')
    parser.add_argument('--manifest', default=None,
                         help='Path to write a JSON manifest of every generated axis-tuple + its jobs_folder/eos path '
                              '(default: <jobs_folder_base>/matrix_manifest.json)')
    args = parser.parse_args()

    batch_name = args.batch_name or f"MassComb_fitmatrix_{datetime.date.today():%Y%m%d}"

    entries = build_matrix()
    if 'all' not in args.categories:
        entries = [e for e in entries if e['category'] in args.categories]

    resolved = []
    skipped = []
    for e in entries:
        if not valid_combo(e['dataset_key'], e['theory_source'], e['order'], e['pdf'], e['poi_config']):
            skipped.append(e)
            continue
        resolved.append(e)

    if skipped:
        print(f"Skipping {len(skipped)} invalid axis-tuple(s) (failed matrix_axes.valid_combo):")
        for e in skipped:
            print(f"  {e}")

    # Cross every resolved axis-tuple with both interpolation-polynomial
    # orders: user decision 2026-07-12 to always get both settings' final
    # results, mirrors matrix_axes.POLY_ORDERS as the single source of truth.
    resolved = [dict(e, poly_order=p) for e in resolved for p in POLY_ORDERS]

    print(f"Generating {len(resolved)} axis-tuple(s) (x{len(POLY_ORDERS)} poly orders) "
          f"across categories {sorted(set(e['category'] for e in resolved))}")

    # Every axis-tuple sharing the same ConvinoSetup produces a BYTE-IDENTICAL
    # tarball (setup dir + convino exe + mtpole-ttj + the full, unfiltered
    # theory JSONs -- --PDF/--order/--polyOrder are runtime doFit.py flags in
    # the generated .sh, never baked into the tarball). createCondorJobs.py
    # namespaces its tarball path by --jobs-folder, and this script gives
    # every axis-tuple its own jobs_folder (needed for --eos-output-path
    # collision-avoidance, see module docstring), so without this cache each
    # of the ~30 axis-tuples per dataset independently re-tars the same ~31MB
    # payload -- confirmed on the real 2026-07-13 run (232 axis-tuples across
    # only 7 unique setups drove local disk to several GB before this fix).
    # Fix: build each unique setup's tarball once, then symlink (not copy --
    # zero extra disk) it plus its .meta.json into every later axis-tuple's
    # jobs_folder before calling createCondorJobs.py with --reuse-tarball, so
    # it recognizes the metadata match and skips recreating it. A real
    # hardlink would be preferable (indistinguishable from a plain file to
    # any reader, including Condor's file transfer) but AFS does not support
    # hard links across directories at all, even within the same volume --
    # confirmed directly (os.link raises EXDEV; os.symlink across the same
    # two directories works fine, and os.path.isfile() follows it correctly).
    tarball_cache = {}

    manifest = []
    for e in resolved:
        tag = axis_tag(e)
        setup = DATASET_KEY_TO_SETUP[e['dataset_key']]
        setup_path = os.path.join(REPO_ROOT, 'pyconvino', 'ConvinoSetups', setup)
        eos_output_path = os.path.join(args.eos_output_base, tag)
        jobs_folder = os.path.join(args.jobs_folder_base, tag)
        dofit_extra_args = dofit_extra_args_for(e, args.common_dofit_extra_args)

        tarball_filename = sanitize(os.path.basename(setup_path)) + '_package.tgz'
        tarball_path = os.path.join(jobs_folder, tarball_filename)
        meta_path = tarball_path + '.meta.json'
        cached = tarball_cache.get(setup)
        if cached and not args.dry_run:
            os.makedirs(jobs_folder, exist_ok=True)
            os.symlink(os.path.abspath(cached[0]), tarball_path)
            shutil.copyfile(cached[1], meta_path)

        cmd = [
            sys.executable, os.path.join(REPO_ROOT, 'condorManagement', 'createCondorJobs.py'),
            setup_path, eos_output_path, f"{batch_name}__{tag}", jobs_folder,
            '--only-nominal', '--mode', 'new',
            '--conda-pack-tarball', args.conda_pack_tarball,
            '--nlo-theory-json', args.nlo_theory_json,
            '--stripper-theory-json', args.stripper_theory_json,
            '--blind-salt-file', args.blind_salt_file,
            '--dofit-extra-args', dofit_extra_args,
            '--reuse-tarball',
        ]
        if args.do_impacts:
            cmd.append('--do-impacts')

        manifest.append({
            'axis_tag': tag, 'category': e['category'], 'dataset_key': e['dataset_key'],
            'theory_source': e['theory_source'], 'order': e['order'], 'pdf': e['pdf'],
            'poi_config': e['poi_config'], 'poly_order': e['poly_order'], 'setup_path': setup_path,
            'eos_output_path': eos_output_path, 'jobs_folder': jobs_folder,
            'dofit_extra_args': dofit_extra_args,
        })

        print(f"[{e['category']}] {tag}" + (' (reusing tarball via symlink)' if cached else ''))
        print('  ' + ' '.join(cmd))
        if not args.dry_run:
            subprocess.run(cmd, check=True)
            if setup not in tarball_cache:
                tarball_cache[setup] = (tarball_path, meta_path)

    manifest_path = args.manifest or os.path.join(args.jobs_folder_base, 'matrix_manifest.json')
    if not args.dry_run:
        os.makedirs(os.path.dirname(manifest_path) or '.', exist_ok=True)
        with open(manifest_path, 'w') as f:
            json.dump({'batch_name': batch_name, 'generated': datetime.datetime.now().isoformat(),
                       'entries': manifest}, f, indent=2)
        print(f"Wrote manifest: {manifest_path}")
    else:
        print(f"[dry-run] would write manifest: {manifest_path}")


if __name__ == '__main__':
    main()
