# write a script that calls createCondorJobs.py with all combinations
#!/bin/bash

# python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/Combination_ATLAS813CMS13_noCorr /eos/cms/store/group/cmst3/group/top/sewuchte/MassCombination/condor/ TestBatchName jobs_dir_Combination_ATLAS813CMS13_noCorr --only-nominal

eosOutPath="/eos/cms/store/group/cmst3/group/top/sewuchte/MassCombination/condor/_$(date +%Y%m%d)/"
batchname="MassComb_Convino__$(date +%Y%m%d)"

# create a list of standaline setups
setups=(
"ATLAS13Only"
"ATLAS8Only"
"CMSOnly"
)

# a list of noCorr setups
noCorrSetups=(
"Combination_ATLAS13CMS13_noCorr"
"Combination_ATLAS813CMS13_noCorr"
"Combination_ATLAS813_noCorr"
"Combination_ATLAS8CMS13_noCorr"
)

# a list of corr setups
corrSetupsNoScan=(
"Combination_ATLAS13CMS13_corr"
"Combination_ATLAS813_corr"
"Combination_ATLAS8CMS13_corr"
)

corrSetupsWithScan=(
"Combination_ATLAS813CMS13_corrV2"
)

corrSetupsNoScanTMP=(
# "Combination_ATLAS813CMS13_corrV2_noLineshape"
"Combination_ATLAS813CMS13_corrV2_noRecoil"
)

toTestAtlas=(
"Combination_ATLAS813CMS13_corrExtreme"
"Combination_ATLAS813CMS13_corrMore"
"Combination_ATLAS813CMS13_corrDeltaPhi"
)

# create a jobsdir with datestamp
jobsdir="condorJobs_$(date +%Y%m%d)"

# check if condorJobs directory exists, if not create it
if [ ! -d "$jobsdir" ]; then
    mkdir "$jobsdir"
fi

# and now call the condor jobs
# for setup in "${setups[@]}"; do
#     echo "Creating condor jobs for standalone setup: ${setup}"
#     python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup}  --only-nominal
# done

# for setup in "${noCorrSetups[@]}"; do
#     echo "Creating condor jobs for noCorr setup: ${setup}"
#     python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} --only-nominal
# done

# for setup in "${corrSetupsNoScan[@]}"; do
#     echo "Creating condor jobs for corr setup without scan: ${setup}"
#     python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} --only-nominal
# done

# for setup in "${corrSetupsWithScan[@]}"; do
#     echo "Creating condor jobs for corr setup with scan: ${setup}"
#     python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup}
# done


# # and now call the condor jobs for impacts
# for setup in "${setups[@]}"; do
#     echo "Creating condor jobs for standalone setup: ${setup}"
#     python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup}  --only-nominal --do-impacts
# done

# for setup in "${noCorrSetups[@]}"; do
#     echo "Creating condor jobs for noCorr setup: ${setup}"
#     python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} --only-nominal --do-impacts
# done

# for setup in "${corrSetupsNoScan[@]}"; do
#     echo "Creating condor jobs for corr setup without scan: ${setup}"
#     python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} --only-nominal --do-impacts
# done

# for setup in "${corrSetupsWithScan[@]}"; do
#     echo "Creating condor jobs for corr setup with scan: ${setup}"
#     python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} --only-nominal --do-impacts
# done


# and the tests for TMP
# for setup in "${corrSetupsNoScanTMP[@]}"; do
#     echo "Creating condor jobs for corr setup without scan (TMP): ${setup}"
#     python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} --only-nominal
# done

for setup in "${toTestAtlas[@]}"; do
    echo "Creating condor jobs for corr setup without scan: ${setup}"
    python3 condorManagement/createCondorJobs.py Convino/ConvinoSetups/${setup} ${eosOutPath} ${batchname} ${jobsdir}/${setup} --only-nominal
done