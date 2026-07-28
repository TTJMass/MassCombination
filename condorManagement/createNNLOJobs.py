#!/usr/bin/env python3
"""Create HTCondor jobs to run apptainer ttbarj tests and copy outputs to EOS.

Generates .sh wrappers and .htc submit files compatible with HTC_monitor.py
(it writes a marker <JOBNAME>.completed into the submit directory when finished).

Usage: createNNLOJobs.py --nJobs N --nJobsNNLO M --runtime 3600 --runtimeNNLO 7200 \
                       --nCores 4 --mTop 172.5 --eos-output /eos/.../outdir \
                       --jobs-folder /path/to/jobs_folder
"""
import os
import stat
import argparse
import textwrap

from createCondorJobs import sanitize

HTC_TEMPLATE = """universe = vanilla
executable = {script_name}
batch_name = {jobname}
initialdir = {jobdir}
when_to_transfer_output = ON_EXIT
output = {jobdir}/out/{jobname}.out
error = {jobdir}/err/{jobname}.err
log = {jobdir}/log/{jobname}.log
request_cpus = {request_cpus}
request_memory = 2000MB
+MaxRuntime = {maxruntime}
+AccountingGroup = "group_u_CMST3.all"
queue
"""

SH_TEMPLATE = textwrap.dedent("""#!/bin/bash
set -euo pipefail

JOBNAME={jobname}
ncores={ncores}
TOTAL_TIME={time}
mt={mt}
ECMS=13000

pwd
ls -alh

rsync -av {sif_path} ./

mkdir -p output

{loop_body}

pwd
ls -alh

# copy outputs to central EOS folder (merge contents into single shared folder)
OUTDIR="{eos_outdir}"
mkdir -p "$OUTDIR"
if [ -d "output" ]; then
    # copy contents (including hidden) of output into the shared OUTDIR
    cp -r output/. "$OUTDIR/" || true
fi

# create a marker file so HTC_monitor.py can detect completion
if [ -n "$JOBNAME" ]; then
    touch "{marker_path}" || true
fi
echo "Job finished, outputs copied to $OUTDIR"
""")

def _distribute_times(total_time, perc_map, order):
    """Distribute total_time (seconds) to keys in order according to perc_map (percent values).

    perc_map: dict key->percent (e.g. 47.0 for 47%)
    order: list of keys in desired order
    Returns dict key->seconds (int), sum equals total_time (mod rounding adjusted).
    """
    per_seconds = {}
    for k in order:
        pct = perc_map.get(k, 0.0)
        per_seconds[k] = int(round(total_time * (pct / 100.0)))
    s = sum(per_seconds.values())
    diff = int(total_time) - s
    if diff != 0:
        # add remainder to DU0 if present, else to first
        key = 'DU0' if 'DU0' in per_seconds else order[0]
        per_seconds[key] = per_seconds.get(key, 0) + diff
    return per_seconds


def _build_normal_sequence(sif_path, total_time, job_index, args=None):
    # order used previously
    order = ['B', 'RF', 'UC0', 'UC1', 'VF', 'RVF', 'DU0', 'FR0', 'SU1', 'DU1', 'FR1']
    # percentages as specified by user
    perc_map = {
        'B': 0.03,
        'UC0': 0.03,
        'UC1': 0.03,
        'VF': 0.03,
        'RF': 0.33,
        'DU1': 0.03,
        'FR1': 0.03,
        'SU1': 0.03,
        'FR0': 0.03,
        'DU0': 47.0,
        'RVF': 47.0,
    }
    per_sec = _distribute_times(total_time, perc_map, order)
    parts = []
    # Seed per subpart: offset by job_index*nCores (so jobs don't overlap,
    # assuming each job uses nCores cores with its own per-core offset in the
    # executable) plus idx*(nCores*nJobs)+idx (so subparts don't overlap
    # across jobs either).
    seed_map = {c: args.seed_offset + int(job_index)*args.nCores + idx*(args.nCores*args.nJobs) + idx
                for idx, c in enumerate(order)}
    for con in order:
        total_con = per_sec.get(con, 0)
        if total_con <= 0:
            total_con = 1
        seed_val = seed_map.get(con)
        parts.append(f"apptainer run --bind ./output:/work/output {sif_path} $ncores $ECMS $mt {con} {total_con} {seed_val}")
        parts.append("")
    return "\n".join(parts)


def _build_nnlo_sequence(sif_path, total_time, job_index, args=None):
    # Distribute equally across the two NNLO channels
    cons = ['RRF', 'SU0']
    per_sec = {}
    base = int(round(total_time / len(cons)))
    s = base * len(cons)
    diff = int(total_time) - s
    per_sec[cons[0]] = base + diff
    per_sec[cons[1]] = base
    parts = []
    seed_map = {cons[0]: args.seed_offset + int(job_index)*args.nCores + 100000,
                cons[1]: args.seed_offset + int(job_index)*args.nCores + 1000000}
    for con in cons:
        total_con = per_sec[con]
        if total_con <= 0:
            total_con = 1
        seed_val = seed_map.get(con)
        parts.append(f"apptainer run --bind ./output:/work/output {sif_path} $ncores $ECMS $mt {con} {total_con} {seed_val}")
        parts.append("")
    return "\n".join(parts)


def write_file(path, content, make_executable=False):
    with open(path, 'w') as f:
        f.write(content)
    if make_executable:
        st = os.stat(path)
        os.chmod(path, st.st_mode | stat.S_IEXEC)


def main():
    parser = argparse.ArgumentParser(description='Create HTCondor jobs running apptainer ttbarj tests')
    parser.add_argument('--nJobs', type=int, default=0, help='Number of regular jobs to create')
    parser.add_argument('--nJobsNNLO', type=int, default=0, help='Number of NNLO jobs to create')
    parser.add_argument('--runtime', type=int, default=3600, help='Runtime (seconds) for regular jobs')
    parser.add_argument('--runtimeNNLO', type=int, default=7200, help='Runtime (seconds) for NNLO jobs')
    parser.add_argument('--nCores', type=int, default=4, help='Number of cores per job')
    parser.add_argument('--mTop', type=float, default=172.5, help='Top quark mass to pass to executable')
    parser.add_argument('--eos-output', dest='eos_output', required=True, help='Central EOS output folder path')
    parser.add_argument('--jobs-folder', dest='jobs_folder', default=os.path.join(os.getcwd(), 'condorJobs_NNLO'), help='Folder where job files will be written')
    parser.add_argument('--sif-path', dest='sif_path', default='/afs/cern.ch/work/s/sewuchte/private/MassCombination/NNLORun/ttbarj-cms-atlas-v2.sif', help='Path to apptainer SIF file to transfer')
    parser.add_argument('--seed-offset', type=int, default=0, help='Offset to add to job index for random seeds (to avoid overlap across different jobs)')
    args = parser.parse_args()

    # if jobs_folder exsists already, print a warning and o nothing
    if os.path.exists(args.jobs_folder):
        print(f"Error: jobs folder {args.jobs_folder} already exists. Please remove it or choose a different folder.")
        return

    jobs_folder = os.path.abspath(args.jobs_folder)
    os.makedirs(jobs_folder, exist_ok=True)
    log_dir = os.path.join(jobs_folder, 'log')
    os.makedirs(log_dir, exist_ok=True)
    err_dir = os.path.join(jobs_folder, 'err')
    os.makedirs(err_dir, exist_ok=True)
    out_dir = os.path.join(jobs_folder, 'out')
    os.makedirs(out_dir, exist_ok=True)

    sif_abspath = os.path.abspath(args.sif_path)
    sif_basename = os.path.basename(sif_abspath)

    # create requested number of regular jobs
    for i in range(args.nJobs):
        jobname = f"theory_ttj_job_lonlo_{i:04d}"
        sh_name = os.path.join(jobs_folder, jobname + '.sh')
        htc_name = os.path.join(jobs_folder, jobname + '.htc')

        sh_content = SH_TEMPLATE.format(
            jobname=jobname,
            ncores=args.nCores,
            time=args.runtime,
            mt=args.mTop,
            sif_path=sif_abspath,
            loop_body=_build_normal_sequence(sif_basename, args.runtime, i, args),
            eos_outdir=args.eos_output+"/output/",
            marker_path=os.path.join(jobs_folder, jobname + '.completed')
        )

        write_file(sh_name, sh_content, make_executable=True)

        htc_content = HTC_TEMPLATE.format(
            script_name=os.path.abspath(sh_name),
            jobname=jobname,
            jobdir=jobs_folder,
            request_cpus=args.nCores,
            maxruntime=(int(args.runtime*10.)) # add some buffer to the requested runtime to avoid premature killing of jobs that run slightly longer than expected
        )
        write_file(htc_name, htc_content)

    # create NNLO jobs
    for j in range(args.nJobsNNLO):
        jobname = f"theory_ttj_job_nnlo_{j:03d}"
        sh_name = os.path.join(jobs_folder, jobname + '.sh')
        htc_name = os.path.join(jobs_folder, jobname + '.htc')

        sh_content = SH_TEMPLATE.format(
            jobname=jobname,
            ncores=args.nCores,
            time=args.runtimeNNLO,
            mt=args.mTop,
            sif_path=sif_abspath,
            loop_body=_build_nnlo_sequence(sif_basename, args.runtimeNNLO, j, args),
            eos_outdir=args.eos_output+"/output/",
            marker_path=os.path.join(jobs_folder, jobname + '.completed')
        )

        write_file(sh_name, sh_content, make_executable=True)

        htc_content = HTC_TEMPLATE.format(
            script_name=os.path.abspath(sh_name),
            jobname=jobname,
            jobdir=jobs_folder,
            request_cpus=args.nCores,
            maxruntime=(int(args.runtimeNNLO * 5.))  # add some buffer to the requested runtime to avoid premature killing of jobs that run slightly longer than expected
        )
        write_file(htc_name, htc_content)

    print(f"Wrote {args.nJobs} regular jobs and {args.nJobsNNLO} NNLO jobs to {jobs_folder}")


if __name__ == '__main__':
    main()
