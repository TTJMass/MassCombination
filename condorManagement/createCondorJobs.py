#!/usr/bin/env python3
"""Create HTCondor jobs that vary correlation parameters from -1.0 to 1.0

Creates per-parameter jobs that:
 - package the ConvinoSetup folder, the convino executable and the mtpole-ttj folder
 - modify the extra_correlations.txt locally to set the chosen correlation value
 - run convino with a unique prefix
 - run the doFit.py script on the produced result file
 - copy results to the specified EOS output path

The script writes a shell wrapper and an HTCondor submit (.htc) file for each job.
Designed to be compatible with HTC_monitor.py (it looks for .htc files).
"""
import argparse
import os
import re
import shutil
import stat
import tarfile
import json
from math import isclose

# Local generated/runtime folders inside mtpole-ttj-pyconvino that Condor workers
# never need (they write their own output on the worker) and that can grow to
# dominate the tarball size (output/output_pre-preapp accumulate every local
# doFit.py run; plots_theory* accumulate every local theory-plot run).
MTP_TARBALL_EXCLUDES = ['plots_fit', 'plots_interp', 'rhoPlotNNLO', 'logs',
                         'output', 'output_pre-preapp', 'plots_theory', '__pycache__']


def parse_extra_correlations(extra_file_path):
    """Parse extra_correlations.txt into a dict of mappings.

    Returns dict where keys are unique ids 'lhs__rhs__i' and values are dicts:
      {'lhs': <lhs>, 'rhs': <rhs>, 'value': <float>, 'line': <original_line>}
    """
    entries = {}
    if not os.path.isfile(extra_file_path):
        raise FileNotFoundError(extra_file_path)

    with open(extra_file_path, 'r') as f:
        lines = f.readlines()

    idx = 0
    pattern = re.compile(r"^\s*([^#\s=]+)\s*=\s*\(([^)]+)\)\s*([^#\s]+)")

    for ln in lines:
        m = pattern.match(ln)
        if m:
            lhs = m.group(1).strip()
            val = float(m.group(2).strip())
            rhs = m.group(3).strip()
            key = f"{lhs}__{rhs}__{idx}"
            entries[key] = {'lhs': lhs, 'rhs': rhs, 'value': val, 'line': ln.rstrip('\n')}
            idx += 1

    return entries


def frange(start, stop, step):
    vals = []
    v = start
    # Avoid floating point drift by using integer steps
    nsteps = int(round((stop - start) / step))
    for i in range(nsteps + 1):
        vals.append(round(start + i * step, 6))
    return vals


def make_tarball(components, tarball_path, exclude_for_dirs=None):
    """Create a tar.gz archive containing the listed files/dirs.

    components: list of absolute paths to files or directories, or (path, arcname)
                tuples to override the default top-level arcname (e.g. to nest a
                file under a subdirectory instead of dropping it at the tarball root)
    exclude_for_dirs: dict mapping component absolute path -> list of substrings;
                      any file/dir whose relative path contains one of those
                      substrings will be skipped when adding that component.
    """
    exclude_for_dirs = exclude_for_dirs or {}

    with tarfile.open(tarball_path, 'w:gz') as tar:
        for comp in components:
            arcname_override = None
            if isinstance(comp, tuple):
                comp, arcname_override = comp
            comp = os.path.abspath(comp)
            base = arcname_override or os.path.basename(comp.rstrip('/'))

            if os.path.isfile(comp):
                tar.add(comp, arcname=base)
                continue

            # directory: walk and add files unless excluded
            excl = exclude_for_dirs.get(comp, [])

            for root, dirs, files in os.walk(comp):
                # compute relative path to comp
                relroot = os.path.relpath(root, comp)
                # skip directories whose rel path contains any excluded substring
                skip_dir = False
                if relroot == '.':
                    relroot = ''
                for s in excl:
                    if s in relroot:
                        skip_dir = True
                        break
                if skip_dir:
                    # prevent walking into these dirs
                    dirs[:] = []
                    continue

                # add directory entry (so empty dirs are preserved)
                arcdir = os.path.join(base, relroot) if relroot else base
                tar.add(root, arcname=arcdir, recursive=False)

                for fname in files:
                    relpath = os.path.join(relroot, fname) if relroot else fname
                    # skip files with excluded substrings in their relative path
                    excluded = any(s in relpath for s in excl)
                    if excluded:
                        continue
                    fullpath = os.path.join(root, fname)
                    arcname = os.path.join(base, relpath)
                    tar.add(fullpath, arcname=arcname)


def single_exp_npz_components():
    """(path, arcname) tuples for the three single-experiment pyconvino NPZ fits
    (CMS-only, ATLAS8-only, ATLAS13-only). plotConstraintsAndPulls.py's
    --pulls/--covCompare reads these via configs.input_f's
    '../pyconvino/out/<name>.npz' paths, resolved relative to mtpole-ttj-pyconvino/
    (see the CWD-independent fallback in _load_single_exp_objs there). Shipped
    nested under pyconvino/out/ -- a sibling of mtpole-ttj-pyconvino/ in the
    extracted tarball -- so that fallback finds them on the worker node, which has
    no other copy of pyconvino/out/.
    """
    pyconvino_out = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'pyconvino', 'out'))
    components = []
    for name in ('CMSOnly_result.npz', 'ATLAS8Only_result.npz', 'ATLAS13Only_result.npz'):
        npz_path = os.path.join(pyconvino_out, name)
        if os.path.isfile(npz_path):
            components.append((npz_path, os.path.join('pyconvino', 'out', name)))
        else:
            print(f"WARNING: single-exp NPZ not found, --covCompare will be incomplete on Condor: {npz_path}")
    return components


def sanitize(s):
    return re.sub(r'[^0-9A-Za-z_.-]', '_', s)

# +MaxRuntime = 1209600


HTC_TEMPLATE = """universe = vanilla
executable = {script_name}
batch_name = {jobname}
initialdir = {jobdir}
should_transfer_files = YES
transfer_input_files = {transfer_input_files}
when_to_transfer_output = ON_EXIT
output = {jobdir}/logs/{jobname}.out
error = {jobdir}/logs/{jobname}.err
log = {jobdir}/logs/{jobname}.log
request_cpus = 4
request_memory = 12000MB
+MaxRuntime = 86400
+AccountingGroup = "group_u_CMST3.all"
queue
"""
# do 2 weeks instead of nextweek
# +JobFlavour = "nextweek"

SH_TEMPLATE = """#!/usr/bin/env bash
set -euo pipefail

JOBNAME={jobname}
PARAM_KEY={param_key}
PARAM_LHS={lhs}
PARAM_RHS={rhs}
SETUP_DIR={setup_dir}
PREFIX={prefix}
NEWVAL={value}
OUTDIR={eos_outdir}
TARBALL={tarball_name}
EXEC_DIR=$(pwd)
TAR_BASENAME={tarball_basename}

# TARBALL may be a relative path (relative to EXEC_DIR); resolve it to an
# absolute path now, before we cd into WORKDIR below, otherwise the relative
# path no longer resolves to the right location.
if [[ "$TARBALL" != /* ]]; then
    TARBALL="$EXEC_DIR/$TARBALL"
fi

WORKDIR=$(mktemp -d /tmp/${{USER}}_${{JOBNAME}}_XXXX)
echo "Working dir: $WORKDIR"
cd "$WORKDIR"

# if the package was transferred by Condor into the job's initial directory,
# copy it into our temporary workdir so we can extract it there
if [ -f "$EXEC_DIR/$TAR_BASENAME" ]; then
    cp "$EXEC_DIR/$TAR_BASENAME" . || true
    TARBALL="./$TAR_BASENAME"
fi

tar -xzf "$TARBALL"

# Optional second tarball (--shared-package-tarball): mtpole-ttj-pyconvino +
# theory JSONs, when shipped separately from the (small, setup-specific)
# $TARBALL above -- see createCondorJobs.py's --shared-package-tarball help.
# Extracts into the same WORKDIR; SHARED_TAR_BASENAME is empty (this whole
# block is a no-op) when the flag wasn't used.
SHARED_TARBALL={shared_tarball_name}
SHARED_TAR_BASENAME={shared_tarball_basename}
if [ -n "$SHARED_TAR_BASENAME" ]; then
    if [[ "$SHARED_TARBALL" != /* ]]; then
        SHARED_TARBALL="$EXEC_DIR/$SHARED_TARBALL"
    fi
    if [ -f "$EXEC_DIR/$SHARED_TAR_BASENAME" ]; then
        SHARED_TARBALL="$EXEC_DIR/$SHARED_TAR_BASENAME"
    fi
    tar -xzf "$SHARED_TARBALL"
fi

# find the extracted setup dir (we expect the tarball to contain the setup folder named as basename of provided setup)
EXTRACTED_SETUP_DIR="$(find . -maxdepth 2 -type d -name '{setup_basename}' -print -quit)"
if [ -z "$EXTRACTED_SETUP_DIR" ]; then
  # fallback: use provided relative path
  EXTRACTED_SETUP_DIR="{setup_basename}"
fi

EXTRA_FILE="$(find "$EXTRACTED_SETUP_DIR" -type f -name extra_correlations.txt -print -quit)"
if [ -z "$EXTRA_FILE" ]; then
  echo "ERROR: extra_correlations.txt not found in setup" >&2
  exit 2
fi

echo "Modifying $EXTRA_FILE: setting $PARAM_LHS -> $NEWVAL"
# only modify the extra_correlations file if NEWVAL is not the special token NOCHANGE
if [ "${{NEWVAL}}" != "NOCHANGE" ]; then
# pass filename, lhs, rhs and new value as argv so shell variables are expanded
python3 - "$EXTRA_FILE" "$PARAM_LHS" "$PARAM_RHS" {value} <<'PY'
import re, sys
fn = sys.argv[1]
lhs = sys.argv[2]
rhs = sys.argv[3]
newv = float(sys.argv[4])
txt = open(fn).read()
pat = re.compile(r"(^\s*"+re.escape(lhs)+r"\s*=\s*\()\s*[^)]+(\)\s*"+re.escape(rhs)+r")", flags=re.M)
def repl(m):
    return m.group(1)+str(newv)+m.group(2)
txt2, n = pat.subn(repl, txt)
if n == 0:
    # try more relaxed match (any rhs)
    pat2 = re.compile(r"(^\s*"+re.escape(lhs)+r"\s*=\s*\()\s*[^)]+(\)\s*.*$)", flags=re.M)
    txt2, n2 = pat2.subn(repl, txt)
    if n2 == 0:
        print('WARNING: did not find explicit line for', lhs)
open(fn,'w').write(txt2)
PY
fi

# Use a self-contained, locally-extracted copy of the masscomb conda env
# instead of activating it from EOS. Under heavy concurrent batch load, EOS
# returns transient "Input/output error" on files read during Python
# interpreter startup (site.py, .pth processing) when thousands of jobs hit
# the same EOS-hosted env at once; this doesn't reproduce locally (single
# job, no contention) but caused most scan jobs to fail on the cluster.
# The packed env is shipped like the setup tarball via Condor's own file
# transfer (from the submit host, not read concurrently by workers off EOS),
# so extracting/activating it here has no EOS runtime dependency at all.
ENV_TARBALL={env_tarball_name}
ENV_TAR_BASENAME={env_tarball_basename}
if [[ "$ENV_TARBALL" != /* ]]; then
    ENV_TARBALL="$EXEC_DIR/$ENV_TARBALL"
fi
if [ -f "$EXEC_DIR/$ENV_TAR_BASENAME" ]; then
    ENV_TARBALL="$EXEC_DIR/$ENV_TAR_BASENAME"
fi

# Blinding salt (see mtpole-ttj-pyconvino/blinding.py): transferred as a plain
# (non-tarball) file, same as the setup/env tarballs above. A dotfile name
# (BLIND_SALT_NAME starts with '.') so it is never swept up by the output-copy
# glob near the end of this script, which doesn't match dotfiles. Empty
# BLIND_SALT_NAME (--mode old, or no --blind-salt-file configured) makes the
# -f check below false and this becomes a no-op.
BLIND_SALT_NAME={blind_salt_basename}
if [ -n "$BLIND_SALT_NAME" ] && [ -f "$EXEC_DIR/$BLIND_SALT_NAME" ]; then
    cp "$EXEC_DIR/$BLIND_SALT_NAME" .
    export MASSCOMB_BLIND_SALT_PATH="$WORKDIR/$BLIND_SALT_NAME"
fi

mkdir -p masscomb_env
tar -xzf "$ENV_TARBALL" -C masscomb_env

# conda-pack's activate script isn't written for `set -euo pipefail`: it does
# `type deactivate; if [ $? -eq 0 ]; then ...` (errexit fires on the `type`
# check itself, before the `if` inspects $?), and it references $PS1 (unset
# under nounset in a non-interactive script). Relax both for this step only.
set +eu
source masscomb_env/bin/activate
conda-unpack
set -euo pipefail

# Never fall back to $HOME/.local/lib/.../site-packages: on worker nodes that
# path lives on AFS, which is frequently unavailable/unreliable in batch jobs.
# All Python dependencies must live inside the extracted masscomb env itself.
export PYTHONNOUSERSITE=1

echo "Running convino"
{convino_cmd} &> convino_run.log || true

# find result file (try common patterns)
RESULT="$(ls ${{PREFIX}}*result*.{result_ext} 2>/dev/null | head -n1 || true)"
if [ -z "$RESULT" ]; then
    RESULT="$(ls *${{PREFIX}}*.{result_ext} 2>/dev/null | grep -i result | head -n1 || true)"
fi

FIT_OK=0
if [ -z "$RESULT" ]; then
  echo "ERROR: No result file found for prefix $PREFIX — convino likely failed; check convino_run.log" >&2
else
    echo "Found result: $RESULT"
    echo "Running doFit.py on $RESULT"
    if python3 {dofit_script}/doFit.py --exp "{exp}" --expPath "$RESULT" --thinputpath "{thinputpath}" {stripperpath_flag} {dofit_extra_args} &> fit.log; then
        FIT_OK=1
        echo "doFit.py succeeded"
    else
        echo "ERROR: doFit.py failed (exit $?) — check fit.log" >&2
    fi
fi

mkdir -p "$OUTDIR"
# always copy logs and any outputs to aid debugging
# "output" (new/pyconvino doFit.py: results.json, budget/interp plots, ...)
# and "plots_fit" (old mtpole-ttj doFit.py) hold the structured fit outputs,
# written to a path relative to $WORKDIR -- without them only the raw stdout
# logs get shipped and results.json / plots never reach EOS.
cp -r convino_run.log {prefix}* *.log *.txt *.pdf plot* plots rhoPlotNNLO output plots_fit "$OUTDIR" 2>/dev/null || true
# only mark completed when the fit actually succeeded
if [ "$FIT_OK" -eq 1 ] && [ -n "$EXEC_DIR" ]; then
    touch "$EXEC_DIR/${{JOBNAME}}.completed" || true
else
    echo "Fit did not succeed — .completed marker NOT written" >&2
fi

echo "Job finished, outputs copied to $OUTDIR"
"""


def main():
    parser = argparse.ArgumentParser(description='Create condor HTCondor jobs varying correlation parameters')
    parser.add_argument('setup_path', help='Path to ConvinoSetup directory (folder containing rho_config.txt and extra_correlations.txt)')
    parser.add_argument('eos_output_path', help='Base output path on EOS (directory)')
    parser.add_argument('batch_name', help='Batch name for jobs')
    parser.add_argument('jobs_folder', help='Folder where job .htc, .sh and logs will be written')
    parser.add_argument('--convino-exe', default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Convino', 'convino')), help='Path to convino executable (default: repository Convino/convino)')
    parser.add_argument('--dofit-path', default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'mtpole-ttj')), help='Path to mtpole-ttj folder containing doFit.py')
    parser.add_argument('--min', type=float, default=-1.0, help='Minimum correlation value')
    parser.add_argument('--max', type=float, default=1.0, help='Maximum correlation value')
    parser.add_argument('--step', type=float, default=0.05, help='Step size')
    parser.add_argument('--max-jobs', type=int, default=None, help='Optional: maximum number of job files to generate (for testing)')
    parser.add_argument('--reuse-tarball', action='store_true', help='If set, reuse existing tarball when components are unchanged')
    parser.add_argument('--force-tarball', action='store_true', help='Force recreation of the tarball even if metadata matches')
    parser.add_argument('--only-nominal', action='store_true', help='Only create the single nominal job and skip creating the scan jobs')
    parser.add_argument('--do-impacts', action='store_true', help='Run also impacts')
    parser.add_argument('--mode', choices=['old', 'new'], default='new',
                        help='Which stack to use: old (Convino C++ binary + mtpole-ttj, decommissioned 2026-07-28 -- see archive/) or new (pyconvino + mtpole-ttj-pyconvino)')
    parser.add_argument('--convino-extra-args', default='',
                        help='Extra arguments appended verbatim to the convino/pyconvino CLI invocation in the generated shell script.')
    parser.add_argument('--dofit-extra-args', default='',
                        help='Extra arguments appended verbatim to the doFit.py invocation in the generated shell script.')
    parser.add_argument('--conda-pack-tarball',
                        default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'envCache', 'masscomb_packed.tar.gz')),
                        help='Path to a conda-pack tarball of the masscomb env, shipped to each job via Condor file '
                             'transfer and extracted locally (avoids activating conda directly from EOS, which is '
                             'unreliable under heavy concurrent batch load). Build it with: '
                             'conda-pack -n masscomb -o <path> (requires no editable-installed packages in the env).')
    parser.add_argument('--nlo-theory-json',
                        default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'theory-data', 'nlo_converted.json')),
                        help='Path to the converted NLO theory JSON (built once offline by '
                             'mtpole-ttj-pyconvino/convertNLOTheoryData.py). When present and --mode new, it is '
                             'shipped inside the setup tarball and used as --thinputpath, so the fit never reads '
                             'the AFS-hosted powheg_generations ROOT files live. Ignored for --mode old.')
    parser.add_argument('--no-ship-nlo-json', action='store_true',
                        help='Force the legacy live-AFS ROOT theory path (inputs/theory_path.txt) even if '
                             '--nlo-theory-json exists. For debugging/fallback only.')
    parser.add_argument('--stripper-theory-json',
                        default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'theory-data', 'withVVF', 'data.json')),
                        help='Path to the converted NNLO/--stripper theory JSON (theory-data/withVVF/data.json). When '
                             'present and --mode new, it is shipped inside the setup tarball and passed as '
                             '--stripperPath, so --stripper fits never read this 115MB file live off AFS '
                             '(doFit.py otherwise defaults --stripperPath to a hardcoded AFS path). Ignored for '
                             '--mode old.')
    parser.add_argument('--no-ship-stripper-json', action='store_true',
                        help='Force the live-AFS --stripperPath default even if --stripper-theory-json exists. '
                             'For debugging/fallback only.')
    parser.add_argument('--blind-salt-file',
                        default=os.path.expanduser('~/.masscomb_blind_salt'),
                        help='Path to the persistent blinding salt file (mtpole-ttj-pyconvino/blinding.py, '
                             'generated once via generate_salt()). Shipped to each job via Condor file '
                             'transfer and exported as MASSCOMB_BLIND_SALT_PATH so blinded (combination) '
                             'fits can find it on the worker. Only relevant for --mode new; ignored for '
                             '--mode old (blinding is not implemented for the legacy stack).')
    parser.add_argument('--shared-package-tarball', default=None,
                        help='Path to a pre-built tarball (see --build-shared-tarball) containing --dofit-path '
                             '(mtpole-ttj-pyconvino) + the theory JSONs -- the part of the per-setup tarball that '
                             'is BYTE-IDENTICAL across every ConvinoSetup. Without this flag, that content is '
                             'bundled into every per-setup tarball as before (backward compatible, just slower '
                             'across multiple distinct setups: each one re-tars the same ~100+MB payload). With '
                             'it, the setup tarball shrinks to just the (small, genuinely setup-specific) '
                             'ConvinoSetup dir, and this shared tarball is shipped/extracted alongside it -- build '
                             'once with --build-shared-tarball, reuse across every setup and every axis-tuple.')
    parser.add_argument('--build-shared-tarball', metavar='OUTPUT_PATH', default=None,
                        help='Build the --shared-package-tarball artifact at OUTPUT_PATH (from --dofit-path + '
                             '--nlo-theory-json + --stripper-theory-json, same exclude rules as the normal setup '
                             'tarball) and exit immediately -- does not generate any jobs. Run this once whenever '
                             'mtpole-ttj-pyconvino source or the theory JSONs change, same lifecycle as '
                             'rebuild_env_pack.sh for the conda env tarball.')
    args = parser.parse_args()

    if args.build_shared_tarball:
        mtp = os.path.abspath(args.dofit_path)
        if args.mode == 'new' and os.path.abspath(mtp) == os.path.abspath(
                os.path.join(os.path.dirname(__file__), '..', 'mtpole-ttj')):
            mtp = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'mtpole-ttj-pyconvino'))
        components = []
        if os.path.isdir(mtp):
            components.append(mtp)
        for json_arg in (args.nlo_theory_json, args.stripper_theory_json):
            json_path = os.path.abspath(json_arg)
            if os.path.isfile(json_path):
                components.append(json_path)
        if args.mode == 'new':
            components.extend(single_exp_npz_components())
        exclude_map = {os.path.abspath(mtp): MTP_TARBALL_EXCLUDES} if os.path.isdir(mtp) else {}
        out_path = os.path.abspath(args.build_shared_tarball)
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        print(f"Building shared package tarball at {out_path} from: {components}")
        make_tarball(components, out_path, exclude_for_dirs=exclude_map)
        print(f"Done: {out_path} ({os.path.getsize(out_path) / 1e6:.1f} MB)")
        return

    conda_pack_tarball = os.path.abspath(args.conda_pack_tarball)
    if not os.path.isfile(conda_pack_tarball):
        raise SystemExit(
            f'--conda-pack-tarball not found: {conda_pack_tarball}\n'
            'Build it once with: conda-pack -n masscomb -o ' + conda_pack_tarball)

    if args.shared_package_tarball:
        args.shared_package_tarball = os.path.abspath(args.shared_package_tarball)
        if not os.path.isfile(args.shared_package_tarball):
            raise SystemExit(
                f'--shared-package-tarball not found: {args.shared_package_tarball}\n'
                'Build it once with: python3 ' + os.path.abspath(__file__) +
                ' --build-shared-tarball ' + args.shared_package_tarball)

    # Set mode-dependent defaults (only if user didn't explicitly pass them)
    if args.mode == 'new':
        old_default = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'mtpole-ttj'))
        if os.path.abspath(args.dofit_path) == old_default:
            args.dofit_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'mtpole-ttj-pyconvino'))

    # do impacts only works with only nominal
    # if args.do_impacts and not args.only_nominal:
    #     raise SystemExit('Error: --do-impacts requires --only-nominal to be set')

    # if args.do_impacts:
    #     # make sure the last "/" is not present
    #     if args.jobs_folder.endswith('/'):
    #         args.jobs_folder = args.jobs_folder[:-1]
    #     if args.eos_output_path.endswith('/'):
    #         print (f"Stripping trailing / from eos_output_path: {args.eos_output_path}")
    #         args.eos_output_path = args.eos_output_path[:-1]
    #     # append _impacts to jobs folder and eos output path
    #     args.jobs_folder = os.path.join(args.jobs_folder, '_impacts')
    #     args.eos_output_path = os.path.join(args.eos_output_path, '_impacts')
    #     print (f"Impacts requested: modified jobs_folder: {args.jobs_folder}, eos_output_path: {args.eos_output_path}")

    setup_path = os.path.abspath(args.setup_path)
    if not os.path.isdir(setup_path):
        raise SystemExit('setup_path not found: ' + setup_path)

    extra_file = os.path.join(setup_path, 'extra_correlations.txt')
    entries = parse_extra_correlations(extra_file)
    if not args.only_nominal and len(entries) == 1:
        if len(entries) == 0:
            raise SystemExit('No correlation mappings found in ' + extra_file)

    os.makedirs(args.jobs_folder, exist_ok=True)
    os.makedirs(os.path.join(args.jobs_folder, 'logs'), exist_ok=True)

    # Compute mode-dependent template values
    dofit_folder = os.path.basename(os.path.abspath(args.dofit_path))

    if args.mode == 'old':
        convino_cmd = './convino -d "$EXTRACTED_SETUP_DIR"/rho_config.txt --prefix "$PREFIX" --noImpacts --neyman'
        if args.do_impacts:
            convino_cmd = './convino -d "$EXTRACTED_SETUP_DIR"/rho_config.txt --prefix "$PREFIX" --neyman'
        if args.convino_extra_args:
            convino_cmd += ' ' + args.convino_extra_args
        result_ext = 'txt'
    else:  # new
        convino_cmd = 'convino "$EXTRACTED_SETUP_DIR"/rho_config.txt --prefix "$PREFIX" --no-impacts --export npz'
        if args.do_impacts:
            convino_cmd = 'convino "$EXTRACTED_SETUP_DIR"/rho_config.txt --prefix "$PREFIX" --export npz'
        if args.convino_extra_args:
            convino_cmd += ' ' + args.convino_extra_args
        result_ext = 'npz'

    # create a tarball of the setup folder + convino exe + mtpole-ttj
    tar_components = [setup_path]
    convino_exe = os.path.abspath(args.convino_exe)
    if args.mode == 'old' and os.path.exists(convino_exe):
        tar_components.append(convino_exe)
    # If --shared-package-tarball is given, mtpole-ttj-pyconvino + the theory
    # JSONs (see below) travel in that separately-shipped, separately-cached
    # tarball instead -- they are BYTE-IDENTICAL across every ConvinoSetup, so
    # baking them into this per-setup tarball means every distinct setup
    # re-tars the same ~100+MB payload for no reason (confirmed: ~10-12min
    # per distinct setup, dominated by AFS I/O in tarfile's directory walk).
    mtp = os.path.abspath(args.dofit_path)
    if os.path.isdir(mtp) and not args.shared_package_tarball:
        tar_components.append(mtp)
        if args.mode == 'new':
            tar_components.extend(single_exp_npz_components())

    # Ship a theory JSON inside the tarball instead of reading it live off AFS
    # from inside the fit (the dominant AFS load source under heavy
    # concurrent scan-job load, see PLAN_afs_load_fix.md). Only the new
    # (pyconvino) stack knows how to read either JSON. Shared by the NLO
    # (theory-data/powheg_generations) and NNLO/--stripper (theory-data/withVVF/data.json,
    # 115MB -- well within what transfer_input_files already proves out at
    # 268MB for the conda-pack env) theory sources below.
    def _ship_or_fallback(json_path, no_ship_flag, fallback_value, flag_name, extra_warning):
        json_path = os.path.abspath(json_path)
        ship = (args.mode == 'new' and not no_ship_flag and os.path.isfile(json_path))
        if args.mode == 'new' and not ship and not no_ship_flag:
            print(f"WARNING: {flag_name} not found at {json_path}; {extra_warning}")
        if ship:
            # Already present in --shared-package-tarball -- just report the
            # basename doFit.py will find after both tarballs are extracted
            # into the same worker workdir, don't duplicate it here too.
            if not args.shared_package_tarball:
                tar_components.append(json_path)
            return os.path.basename(json_path)
        return fallback_value

    thinputpath_value = _ship_or_fallback(
        args.nlo_theory_json, args.no_ship_nlo_json,
        fallback_value=os.path.join(dofit_folder, 'inputs', 'theory_path.txt'),
        flag_name='--nlo-theory-json',
        extra_warning='falling back to live-AFS ROOT theory path (inputs/theory_path.txt). '
                       'Run mtpole-ttj-pyconvino/convertNLOTheoryData.py once to avoid AFS load under heavy batch load.')

    stripperpath_value = _ship_or_fallback(
        args.stripper_theory_json, args.no_ship_stripper_json,
        fallback_value=os.path.abspath(args.stripper_theory_json),
        flag_name='--stripper-theory-json',
        extra_warning="--stripper jobs will fall back to doFit.py's hardcoded live-AFS default.")
    # --stripperPath only exists on the new (pyconvino) doFit.py -- the old
    # (mtpole-ttj) doFit.py has no such argument and argparse hard-errors on
    # any unrecognized flag, so this must be omitted entirely for --mode old.
    stripperpath_flag = f'--stripperPath "{stripperpath_value}"' if args.mode == 'new' else ''

    blind_salt_file = os.path.abspath(os.path.expanduser(args.blind_salt_file))
    ship_blind_salt = (args.mode == 'new' and os.path.isfile(blind_salt_file))
    if args.mode == 'new' and not ship_blind_salt:
        print(f"WARNING: --blind-salt-file not found at {blind_salt_file}; any combination "
              "fit in this batch will hard-fail in blinding._load_salt() on the worker "
              "(unless --unblind is passed via --dofit-extra-args). Run "
              "blinding.generate_salt() once to create it.")
    blind_salt_basename = os.path.basename(blind_salt_file) if ship_blind_salt else ''

    tarball_name = os.path.join(args.jobs_folder, sanitize(os.path.basename(setup_path)) + '_package.tgz')
    # exclude heavy or runtime/generated folders from mtpole-ttj to keep package small
    exclude_map = {}
    if os.path.isdir(mtp):
        exclude_map[os.path.abspath(mtp)] = MTP_TARBALL_EXCLUDES

    # decide whether to (re)create the tarball. If requested, reuse an existing
    # tarball when the set of components and exclude-map keys match previous run.
    meta_path = tarball_name + '.meta.json'
    current_meta = {
        'components': [os.path.abspath(p) for p in tar_components],
        'exclude_keys': sorted(list(exclude_map.keys())),
        'convino_exe': convino_exe,
        'mtp': mtp,
    }

    recreate = True
    if args.reuse_tarball and os.path.isfile(tarball_name) and os.path.isfile(meta_path) and not args.force_tarball:
        try:
            with open(meta_path, 'r') as mf:
                prev = json.load(mf)
            # compare relevant fields
            if prev.get('components') == current_meta['components'] and prev.get('exclude_keys') == current_meta['exclude_keys'] and prev.get('convino_exe') == current_meta['convino_exe'] and prev.get('mtp') == current_meta['mtp']:
                recreate = False
                print(f"Reusing existing tarball {tarball_name} (components unchanged)")
        except Exception:
            recreate = True

    if recreate:
        print(f"Creating tarball {tarball_name} ...")
        make_tarball(tar_components, tarball_name, exclude_for_dirs=exclude_map)
        # write metadata
        try:
            with open(meta_path, 'w') as mf:
                json.dump(current_meta, mf, indent=2)
        except Exception as e:
            print('Warning: could not write tarball metadata:', e)

    values = frange(args.min, args.max, args.step)

    job_count = 0

    if not args.only_nominal:
        for key, info in entries.items():
            lhs = info['lhs']
            rhs = info['rhs']
            # create a correlation-specific subfolder inside the jobs folder
            parent_jobdir = os.path.abspath(args.jobs_folder)
            corr_name = sanitize(f"{lhs}__{rhs}")
            corr_dir = os.path.join(parent_jobdir, corr_name)
            os.makedirs(corr_dir, exist_ok=True)
            os.makedirs(os.path.join(corr_dir, 'logs'), exist_ok=True)

            for val in values:
                jobname = sanitize(f"{os.path.basename(setup_path)}__{lhs}__{rhs}__{val}")
                prefix = jobname
                # create one subfolder per job inside the correlation folder
                jobdir = os.path.join(corr_dir, jobname)
                os.makedirs(jobdir, exist_ok=True)
                os.makedirs(os.path.join(jobdir, 'logs'), exist_ok=True)
                sh_path = os.path.join(jobdir, jobname + '.sh')
                htc_path = os.path.join(jobdir, jobname + '.htc')

                # compute relative tarball paths from this jobdir so the wrapper can find them
                rel_tar = os.path.relpath(tarball_name, start=jobdir)
                rel_env_tar = os.path.relpath(conda_pack_tarball, start=jobdir)
                rel_shared_tar = os.path.relpath(args.shared_package_tarball, start=jobdir) if args.shared_package_tarball else ''
                shared_tarball_basename = os.path.basename(args.shared_package_tarball) if args.shared_package_tarball else ''

                with open(sh_path, 'w') as shf:
                    shf.write(SH_TEMPLATE.format(jobname=jobname,
                                                 param_key=key,
                                                 lhs=lhs,
                                                 rhs=rhs,
                                                 setup_dir=setup_path,
                                                 eos_outdir=os.path.join(os.path.abspath(args.eos_output_path), os.path.basename(setup_path), sanitize(f"{lhs}__{rhs}"), str(val)),
                                                tarball_name=rel_tar,
                                                tarball_basename=os.path.basename(tarball_name),
                                                env_tarball_name=rel_env_tar,
                                                env_tarball_basename=os.path.basename(conda_pack_tarball),
                                                shared_tarball_name=rel_shared_tar,
                                                shared_tarball_basename=shared_tarball_basename,
                                                 setup_basename=os.path.basename(setup_path),
                                                 value=repr(val),
                                                 prefix=prefix,
                                                 exp=sanitize(f"{args.batch_name}__{jobname}"),
                                                 thinputpath=thinputpath_value,
                                                 stripperpath_flag=stripperpath_flag,
                                                 convino_cmd=convino_cmd,
                                                 result_ext=result_ext,
                                                 dofit_script=dofit_folder,
                                                 dofit_extra_args=args.dofit_extra_args,
                                                 blind_salt_basename=blind_salt_basename))

                # make executable
                st = os.stat(sh_path)
                os.chmod(sh_path, st.st_mode | stat.S_IEXEC)

                # HTCondor submit file
                # use absolute paths so condor_submit can be called from any cwd
                abs_tar = os.path.abspath(tarball_name)
                abs_sh = os.path.abspath(sh_path)
                transfer_inputs = abs_tar + ',' + conda_pack_tarball
                if args.shared_package_tarball:
                    transfer_inputs += ',' + args.shared_package_tarball
                if ship_blind_salt:
                    transfer_inputs += ',' + blind_salt_file

                with open(htc_path, 'w') as htf:
                    htf.write(HTC_TEMPLATE.format(script_name=abs_sh,
                                                  transfer_input_files=transfer_inputs,
                                                  jobdir=jobdir,
                                                  jobname=jobname))

                # copy tarball and shell to job folder (they will be transfer_input_files)
                # do not duplicate the tarball per-correlation; condor will transfer the
                # package referenced relatively. Only ensure the shell script is present
                dst_sh = os.path.join(jobdir, os.path.basename(sh_path))
                if os.path.abspath(sh_path) != os.path.abspath(dst_sh):
                    shutil.copy(sh_path, dst_sh)

                job_count += 1
                if args.max_jobs and job_count >= args.max_jobs:
                    print(f"Reached max-jobs={args.max_jobs}; stopping generation")
                    return

    # Create a single global nominal job that uses the input correlations as-is
    # (do not modify extra_correlations.txt). We set NEWVAL to the sentinel
    # NOCHANGE so the wrapper skips the edit step.
    # if we run only-nominal the entry doesnt exist, so just write anything dummy to lhs0/rhs0
    lhs0 = 'DUMMY_LHS'
    rhs0 = 'DUMMY_RHS'
    # first_entry = next(iter(entries.values()))
    # lhs0 = first_entry['lhs']
    # rhs0 = first_entry['rhs']

    parent_jobdir = os.path.abspath(args.jobs_folder)
    nominal_dir = os.path.join(parent_jobdir, 'nominal')
    os.makedirs(nominal_dir, exist_ok=True)
    os.makedirs(os.path.join(nominal_dir, 'logs'), exist_ok=True)

    jobname = sanitize(f"{os.path.basename(setup_path)}__nominal")
    prefix = jobname
    jobdir = os.path.join(nominal_dir, jobname)
    os.makedirs(jobdir, exist_ok=True)
    os.makedirs(os.path.join(jobdir, 'logs'), exist_ok=True)

    sh_path = os.path.join(jobdir, jobname + '.sh')
    htc_path = os.path.join(jobdir, jobname + '.htc')

    rel_tar = os.path.relpath(tarball_name, start=jobdir)
    rel_env_tar = os.path.relpath(conda_pack_tarball, start=jobdir)
    rel_shared_tar = os.path.relpath(args.shared_package_tarball, start=jobdir) if args.shared_package_tarball else ''
    shared_tarball_basename = os.path.basename(args.shared_package_tarball) if args.shared_package_tarball else ''

    with open(sh_path, 'w') as shf:
        if args.do_impacts:
            print("Creating nominal job with impacts")
        shf.write(SH_TEMPLATE.format(jobname=jobname,
                                    param_key='NOMINAL',
                                    lhs=lhs0,
                                    rhs=rhs0,
                                    setup_dir=setup_path,
                                    eos_outdir=os.path.join(os.path.abspath(args.eos_output_path), os.path.basename(setup_path), 'nominal'),
                                    tarball_name=rel_tar,
                                    tarball_basename=os.path.basename(tarball_name),
                                    env_tarball_name=rel_env_tar,
                                    env_tarball_basename=os.path.basename(conda_pack_tarball),
                                    shared_tarball_name=rel_shared_tar,
                                    shared_tarball_basename=shared_tarball_basename,
                                    setup_basename=os.path.basename(setup_path),
                                    value='NOCHANGE',
                                    prefix=prefix,
                                    exp=sanitize(f"{args.batch_name}__{jobname}"),
                                    thinputpath=thinputpath_value,
                                    stripperpath_flag=stripperpath_flag,
                                    convino_cmd=convino_cmd,
                                    result_ext=result_ext,
                                    dofit_script=dofit_folder,
                                    dofit_extra_args=args.dofit_extra_args,
                                    blind_salt_basename=blind_salt_basename))


    st = os.stat(sh_path)
    os.chmod(sh_path, st.st_mode | stat.S_IEXEC)

    abs_tar = os.path.abspath(tarball_name)
    abs_sh = os.path.abspath(sh_path)
    transfer_inputs = abs_tar + ',' + conda_pack_tarball
    if args.shared_package_tarball:
        transfer_inputs += ',' + args.shared_package_tarball
    if ship_blind_salt:
        transfer_inputs += ',' + blind_salt_file
    with open(htc_path, 'w') as htf:
        htf.write(HTC_TEMPLATE.format(script_name=abs_sh,
                                transfer_input_files=transfer_inputs,
                                jobdir=jobdir,
                                jobname=jobname))

    dst_sh = os.path.join(jobdir, os.path.basename(sh_path))
    if os.path.abspath(sh_path) != os.path.abspath(dst_sh):
        shutil.copy(sh_path, dst_sh)

    job_count += 1

    print(f"Written {job_count} job .htc/.sh pairs into {args.jobs_folder} (tarball: {tarball_name})")


if __name__ == '__main__':
    main()
