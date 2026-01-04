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

    components: list of absolute paths to files or directories
    exclude_for_dirs: dict mapping component absolute path -> list of substrings;
                      any file/dir whose relative path contains one of those
                      substrings will be skipped when adding that component.
    """
    exclude_for_dirs = exclude_for_dirs or {}

    with tarfile.open(tarball_path, 'w:gz') as tar:
        for comp in components:
            comp = os.path.abspath(comp)
            base = os.path.basename(comp.rstrip('/'))

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


def sanitize(s):
    return re.sub(r'[^0-9A-Za-z_.-]', '_', s)


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
request_cpus = 1
request_memory = 2000MB
+JobFlavour = "nextweek"
+AccountingGroup = "group_u_CMST3.all"
queue
"""


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

echo "Running convino"
./convino -d "$EXTRACTED_SETUP_DIR"/rho_config.txt --prefix "$PREFIX" --neyman &> convino_run.log || true

# find result file (try common patterns)
RESULT="$(ls ${{PREFIX}}*result*.txt 2>/dev/null | head -n1 || true)"
if [ -z "$RESULT" ]; then
    RESULT="$(ls *${{PREFIX}}*.txt 2>/dev/null | grep -i result | head -n1 || true)"
fi
if [ -z "$RESULT" ]; then
  echo "No result file found for prefix $PREFIX" >&2
else
    echo "Found result: $RESULT"
    echo "Running doFit.py on $RESULT"
    # initialize conda environment using the repository-local conda installation
    if [ -f "/eos/home-s/sewuchte/MyConda/etc/profile.d/conda.sh" ]; then
        source "/eos/home-s/sewuchte/MyConda/etc/profile.d/conda.sh" || true
    else
        export PATH="/eos/home-s/sewuchte/MyConda/bin:$PATH"
    fi
    conda activate masscomb || true
    python3 mtpole-ttj/doFit.py --exp "{exp}" --expPath "$RESULT" --thinputpath "{thinputpath}" &> fit.log || true
fi

mkdir -p "$OUTDIR"
# copy logs, results, pdfs and plot folders
cp -r convino_run.log {prefix}* *.log *.txt *.pdf plot* plots rhoPlotNNLO "$OUTDIR" 2>/dev/null || true
# create a marker file so HTC_monitor.py can detect completion
if [ -n "$EXEC_DIR" ]; then
    touch "$EXEC_DIR/${{JOBNAME}}.completed" || true
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
    args = parser.parse_args()

    setup_path = os.path.abspath(args.setup_path)
    if not os.path.isdir(setup_path):
        raise SystemExit('setup_path not found: ' + setup_path)

    extra_file = os.path.join(setup_path, 'extra_correlations.txt')
    entries = parse_extra_correlations(extra_file)
    if len(entries) == 0:
        raise SystemExit('No correlation mappings found in ' + extra_file)

    os.makedirs(args.jobs_folder, exist_ok=True)
    os.makedirs(os.path.join(args.jobs_folder, 'logs'), exist_ok=True)

    # create a tarball of the setup folder + convino exe + mtpole-ttj
    tar_components = [setup_path]
    convino_exe = os.path.abspath(args.convino_exe)
    if os.path.exists(convino_exe):
        tar_components.append(convino_exe)
    mtp = os.path.abspath(args.dofit_path)
    if os.path.isdir(mtp):
        tar_components.append(mtp)

    tarball_name = os.path.join(args.jobs_folder, sanitize(os.path.basename(setup_path)) + '_package.tgz')
    # exclude heavy or runtime/generated folders from mtpole-ttj to keep package small
    exclude_map = {}
    if os.path.isdir(mtp):
        exclude_map[os.path.abspath(mtp)] = ['plots_fit', 'plots_interp', 'rhoPlotNNLO', 'logs']

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

            # compute relative tarball path from this jobdir so the wrapper can find it
            rel_tar = os.path.relpath(tarball_name, start=jobdir)

            with open(sh_path, 'w') as shf:
                shf.write(SH_TEMPLATE.format(jobname=jobname,
                                             param_key=key,
                                             lhs=lhs,
                                             rhs=rhs,
                                             setup_dir=setup_path,
                                             eos_outdir=os.path.join(os.path.abspath(args.eos_output_path), os.path.basename(setup_path), sanitize(f"{lhs}__{rhs}"), str(val)),
                                            tarball_name=rel_tar,
                                            tarball_basename=os.path.basename(tarball_name),
                                             setup_basename=os.path.basename(setup_path),
                                             value=repr(val),
                                             prefix=prefix,
                                             exp=sanitize(f"{args.batch_name}__{jobname}"),
                                             thinputpath=os.path.join('mtpole-ttj','inputs','theory_path.txt')))

            # make executable
            st = os.stat(sh_path)
            os.chmod(sh_path, st.st_mode | stat.S_IEXEC)

            # HTCondor submit file
            # use absolute paths so condor_submit can be called from any cwd
            abs_tar = os.path.abspath(tarball_name)
            abs_sh = os.path.abspath(sh_path)
            transfer_inputs = abs_tar

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
    first_entry = next(iter(entries.values()))
    lhs0 = first_entry['lhs']
    rhs0 = first_entry['rhs']

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

    with open(sh_path, 'w') as shf:
        shf.write(SH_TEMPLATE.format(jobname=jobname,
                                     param_key='NOMINAL',
                                     lhs=lhs0,
                                     rhs=rhs0,
                                     setup_dir=setup_path,
                                     eos_outdir=os.path.join(os.path.abspath(args.eos_output_path), os.path.basename(setup_path), 'nominal'),
                                     tarball_name=rel_tar,
                                     tarball_basename=os.path.basename(tarball_name),
                                     setup_basename=os.path.basename(setup_path),
                                     value='NOCHANGE',
                                     prefix=prefix,
                                     exp=sanitize(f"{args.batch_name}__{jobname}"),
                                     thinputpath=os.path.join('mtpole-ttj','inputs','theory_path.txt')))

    st = os.stat(sh_path)
    os.chmod(sh_path, st.st_mode | stat.S_IEXEC)

    abs_tar = os.path.abspath(tarball_name)
    abs_sh = os.path.abspath(sh_path)
    transfer_inputs = abs_tar
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
