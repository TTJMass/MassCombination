#!/bin/bash
# Generates HTCondor jobs for the NEW pyconvino + mtpole-ttj-pyconvino stack.
# Mirror of writeAllCondorJobs.sh for the old Convino+mtpole-ttj stack.

eosOutPath="/eos/cms/store/group/cmst3/group/top/sewuchte/MassCombination/condor/_$(date +%Y%m%d)_pyconvino/"
batchname="MassComb_pyconvino__$(date +%Y%m%d)"

setups=(
"ATLAS13Only"
"ATLAS8Only"
"CMSOnly"
)

noCorrSetups=(
"Combination_ATLAS13CMS13_noCorr"
"Combination_ATLAS813CMS13_noCorr"
"Combination_ATLAS813_noCorr"
"Combination_ATLAS8CMS13_noCorr"
)

corrSetupsNoScan=(
"Combination_ATLAS13CMS13_corr"
"Combination_ATLAS813_corr"
"Combination_ATLAS8CMS13_corr"
)

corrSetupsWithScan=(
"Combination_ATLAS813CMS13_corrV2"
)

corrSetupsNoScanTMP=(
"Combination_ATLAS813CMS13_corrV2_noLineshape"
"Combination_ATLAS813CMS13_corrV2_noRecoil"
)

toTestAtlas=(
"Combination_ATLAS813CMS13_corrExtreme"
"Combination_ATLAS813CMS13_corrMore"
"Combination_ATLAS813CMS13_corrDeltaPhi"
)

jobsdir="condorJobs_$(date +%Y%m%d)_pyconvino"

if [ ! -d "$jobsdir" ]; then
    mkdir "$jobsdir"
fi

# Uncomment the blocks you want to run.
# escape properly for the usage of ""
commonArgs=(
    --mode new
    --dofit-extra-args "--interpCheck --nuisanceFit --pulls --budget"
    --do-impacts
    --conda-pack-tarball /afs/cern.ch/work/s/sewuchte/private/MassCombination/envCache/masscomb_packed.tar.gz
)

echo "Creating condor jobs for pyconvino setups..."
echo commonArgs: "${commonArgs[@]}"

for setup in "${setups[@]}"; do
    echo "Creating condor jobs for standalone setup: ${setup}"
    python3 condorManagement/createCondorJobs.py pyconvino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} --only-nominal "${commonArgs[@]}"
done

for setup in "${noCorrSetups[@]}"; do
    echo "Creating condor jobs for noCorr setup: ${setup}"
    python3 condorManagement/createCondorJobs.py pyconvino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} "${commonArgs[@]}" --only-nominal
done

for setup in "${corrSetupsNoScan[@]}"; do
    echo "Creating condor jobs for corr setup without scan: ${setup}"
    python3 condorManagement/createCondorJobs.py pyconvino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} "${commonArgs[@]}" --only-nominal
done

for setup in "${corrSetupsWithScan[@]}"; do
    echo "Creating condor jobs for corr setup with scan: ${setup}"
    python3 condorManagement/createCondorJobs.py pyconvino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} "${commonArgs[@]}"
done

for setup in "${toTestAtlas[@]}"; do
    echo "Creating condor jobs for corr setup without scan: ${setup}"
    python3 condorManagement/createCondorJobs.py pyconvino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} "${commonArgs[@]}" --only-nominal
done
