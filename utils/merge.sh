#!/bin/env bash
set -xeuo pipefail

# rm -rf <dir>

DIRS=( checkperiodic chromatogram extractRNG extractTS genmprd gennmat gmsh2mixd gmsh2mixdv2 mapflow mixd2pvtu mixdclass mixdsplit mixdvolume normmixd rmmat scalemixd shiftmixd stitchperiodic )

# DIRS=( mixdtools )

mkdir -p mixdtools

for DIR in ${DIRS[@]}; do 
    echo "Processing $DIR"
    echo "---------------"
    (git clone git@jugit.fz-juelich.de:IBG-1/ModSim/Chroma-HD/$DIR.git
    cd "$DIR"
    git filter-repo --to-subdirectory-filter $DIR 
    )
    (cd mixdtools
    git remote add $DIR ../$DIR
    git fetch $DIR --no-tags
    EDITOR=true git merge --allow-unrelated-histories $DIR/master
    git remote remove $DIR
    )
    echo "---------------"
done

