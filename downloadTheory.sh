# /bin/sh

function mycopy {
    if [[ $# -lt 2 ]] ; then
        echo "Description: copy files using rsync."  >&2
        echo "Usage: [rsync-options] source destination"  >&2
        exit 1
    fi

    NETWORK_ERRORS=( 12 30 255 )
    while [[ 1 ]] ; do
        rsync --archive --compress --verbose --human-readable --progress --inplace --delete "$@"
        RESULT=$?
        echo "Exit code $RESULT"
        if [[ ! " ${NETWORK_ERRORS[@]} " =~ " $RESULT " ]] ; then break ; fi
        sleep 4
    done
}

mycopy lxplus.cern.ch:/eos/home-s/sewuchte/TopMassCombinationFromTTJ/powheg_generations/ powheg_generations/
mycopy lxplus.cern.ch:/eos/home-s/sewuchte/TopMassCombinationFromTTJ/powheg_generations_8TeV/ powheg_generations_8TeV/