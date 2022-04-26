#!/usr/bin/env bash

set -euo pipefail

## NOTE: Usage: 
# ./clone-all-tools.sh will clone everything listed below to the current directory
#       --https: will pull via https
#       --ssh: will pull via ssh, default
#       any positional arguments left are considered repo names to be pulled

BRANCH="master"
ARGS=""
MODE="SSH"
RMGIT=false # flag to remove the .git folder. Useful in JURECA to not get prompted the password thanks to my shell trying to access git stats
ALL=( checkperiodic extractRNG extractTS genmesh genmprd gennmat gmsh2mixd gmsh2mixdv2 mapflow mixd2pvtu mixdclass mixdsplit mixdvolume normmixd pymesh rmmat scalemixd scripts shiftmixd stitchperiodic xcad )

REPOS=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        --https)
            echo "Cloning with HTTPS"
            MODE="HTTPS"
            shift # past value
            ;;
        --ssh)
            echo "Cloning with SSH"
            MODE="SSH"
            shift # past value
            ;;
        --no-git)
            echo "Performing a shallow clone and removing git history"
            ARGS="$ARGS --depth 1"
            RMGIT=true
            shift # past value
            ;;
        *)    # unknown option
            REPOS+=("$1") # save it in an array for later
            shift # past argument
            ;;
    esac
done
set -- "${REPOS[@]}" # restore positional parameters

if [[ "$MODE" == "SSH" ]]; then 
    SERVER="git@jugit.fz-juelich.de:IBG-1/ModSim/Chroma-HD"
else
    SERVER="https://jugit.fz-juelich.de/IBG-1/ModSim/Chroma-HD"
fi

[[ ${#REPOS[@]} == 0 ]] && REPOS=${ALL[@]}

for REPO in ${REPOS[@]}; do 
    git clone $ARGS "$SERVER/$REPO.git"
    if [ $RMGIT = true ]; then
        rm -rf "$REPO/.git"
    fi
done

# https://jugit.fz-juelich.de/IBG-1/ModSim/Chroma-HD/$REPO/-/archive/master/$REPO-master.tar.gz
