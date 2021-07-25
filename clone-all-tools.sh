#!/usr/bin/env bash

## NOTE: Usage: 
# ./clone-all-tools.sh will clone everything listed below to the current directory
#       --https: will pull via https
#       --ssh: will pull via ssh, default
#       any positional arguments left are considered repo names to be pulled

MODE="SSH"
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
    git clone "$SERVER/$REPO.git"
done
